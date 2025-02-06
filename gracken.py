import os
import sys
import glob
import argparse
import pandas as pd
from ete3 import Tree as Tree3

pd.options.mode.copy_on_write = True


def harmonize_species(name):
    # always remove the "s__" prefix (if present)
    if name.startswith("s__"):
        name = name[3:]
    # replace spaces with underscores if allowed
    # (ape's read.tree function doesn't like spaces in tip names)
    if not args.dont_replace_spaces:
        name = name.replace(" ", "_")
    return name


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a phylogenetic tree using Bracken/Kraken2 reports by pruning the GTDB tree"
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing Bracken/Kraken2 report files",
    )
    parser.add_argument(
        "--bac_taxonomy",
        required=True,
        help="Path to bacterial taxonomy file (e.g., bac120_taxonomy_r95.tsv)",
    )
    parser.add_argument(
        "--ar_taxonomy",
        required=True,
        help="Path to archaeal taxonomy file (e.g., ar122_taxonomy_r95.tsv)",
    )
    parser.add_argument(
        "--bac_tree",
        required=True,
        help="Path to bacterial tree file (e.g., bac120_r95.tree)",
    )
    parser.add_argument(
        "--ar_tree",
        required=True,
        help="Path to archaeal tree file (e.g., ar122_r95.tree)",
    )
    parser.add_argument(
        "--out_prefix",
        default="output",
        help="Prefix for output files (default: output). Creates <prefix>.tree and <prefix>.otu.csv",
    )
    parser.add_argument(
        "--mode",
        choices=["bracken", "kraken2"],
        default="bracken",
        help="Input file format (default: bracken)",
    )
    parser.add_argument(
        "--dont-replace-spaces",
        action="store_true",
        default=False,
        help="Keep spaces in species names (default: False)",
    )
    return parser.parse_args()


def extract_abundances(file_path, mode="bracken"):
    """Extract species abundances from Bracken or Kraken2 reports."""
    df = pd.read_csv(file_path, sep="\t", header=None)

    if mode == "bracken":
        # extract rows where rank is 'S' (species) and abundance > 0
        species_df = df[(df[3] == "S") & (df[2] > 0)]
        species_abundance = species_df[[5, 2]]
    else:
        species_df = df[(df[5] == "S") & (df[1] > 0)]
        species_abundance = species_df[[7, 3]]

    species_abundance.columns = ["species", "abundance"]
    species_abundance["species"] = species_abundance["species"].str.strip()

    return species_abundance


def process_taxonomy(taxonomy_file):
    """Process GTDB taxonomy file and return species mappings."""
    species_to_genomes = {}
    genome_to_species = {}
    with open(taxonomy_file, "r") as f:
        for line in f:
            genome_id, taxonomy = line.strip().split("\t")
            species = taxonomy.split(";")[-1]
            genome_to_species[genome_id] = species
            if species not in species_to_genomes:
                species_to_genomes[species] = []
            species_to_genomes[species].append(genome_id)
    return species_to_genomes, genome_to_species


def find_matching_species(species_set, species_to_genomes):
    """Find matching species."""
    found_species = set()
    missing_species = set()

    for species in species_set:
        if species in species_to_genomes:
            found_species.add(species)
        else:
            missing_species.add(species)

    return found_species, missing_species


def process_tree(tree_file, genome_to_species, species_to_keep):
    """Load, rename and prune tree."""
    tree = Tree3(tree_file, format=1, quoted_node_names=True)

    # replace tip names (genome ids) with species names from taxonomy
    for node in tree.traverse():
        if node.name in genome_to_species:
            node.name = genome_to_species[node.name]

    # harmonize leaf names
    for leaf in tree.get_leaves():
        leaf.name = harmonize_species(leaf.name)

    # get harmonized species names from the tree
    tree_species = set(leaf.name for leaf in tree.get_leaves())
    harmonized_keep = set(harmonize_species(x) for x in species_to_keep)
    valid_species = list(tree_species.intersection(harmonized_keep))

    if not valid_species:
        print(
            f"No species matching between {args.mode} report files and provided taxonomy"
        )
        print("Tree species examples:", list(tree_species)[:5])
        print("Species to keep examples:", list(harmonized_keep)[:5])
        raise ValueError("No matching species found in tree")

    tree.prune(valid_species, preserve_branch_length=True)

    return tree


args = parse_args()


def main():
    # load data and create OTU table
    otu_table = pd.DataFrame()
    file_pattern = "*.breport" if args.mode == "bracken" else "*.k2report"

    for f in glob.glob(os.path.join(args.input_dir, file_pattern)):
        name = os.path.splitext(os.path.basename(f))[0]
        sample_otu = extract_abundances(f, mode=args.mode)

        # Check if sample_otu is empty or doesn't have the required columns
        if sample_otu.empty or not all(
            col in sample_otu.columns for col in ["species", "abundance"]
        ):
            print(
                f"Warning: Skipping {f} because it's empty or missing required columns.",
                file=sys.stderr,
            )
            continue

        sample_otu["sample"] = name
        otu_table = pd.concat([otu_table, sample_otu], ignore_index=True)

    # transform to wide format and harmonize species names
    wide_otu_table = otu_table.pivot(
        index="species", columns="sample", values="abundance"
    )

    # remove 's__' prefix from species names
    wide_otu_table.index = wide_otu_table.index.str.replace(r"^s__", "", regex=True)

    if not args.dont_replace_spaces:
        wide_otu_table.index = wide_otu_table.index.str.replace(" ", "_", regex=False)

    wide_otu_table.to_csv(f"{args.out_prefix}.otu.csv", sep=",")

    # unique species (use original names to match taxonomy)
    species_set = set(otu_table["species"])

    # load taxonomies
    bac_species_to_genomes, bac_genome_to_species = process_taxonomy(args.bac_taxonomy)
    ar_species_to_genomes, ar_genome_to_species = process_taxonomy(args.ar_taxonomy)

    # compare our species against those in the taxonomy
    my_species = {"s__" + x if not x.startswith("s__") else x for x in species_set}
    bac_found, bac_missing = find_matching_species(my_species, bac_species_to_genomes)
    ar_found, ar_missing = find_matching_species(my_species, ar_species_to_genomes)

    # quality control -- make sure we cover all species
    assert ar_missing == bac_found
    assert bac_missing == ar_found

    # make the trees
    bac_tree = process_tree(args.bac_tree, bac_genome_to_species, bac_found)
    ar_tree = process_tree(args.ar_tree, ar_genome_to_species, ar_found)

    # and combine
    combined_tree = Tree3()
    combined_tree.name = "root"
    combined_tree.add_child(bac_tree)
    combined_tree.add_child(ar_tree)

    combined_tree.write(
        outfile=f"{args.out_prefix}.tree", format=1, quoted_node_names=False
    )


if __name__ == "__main__":
    main()
