import os
import sys
import glob
import argparse
import pandas as pd
import gtdb_utils as gt
import ncbi_utils as nc
from ete3 import Tree as Tree3
from data_loaders import read_bracken_report, read_kraken2_report

pd.options.mode.copy_on_write = True


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a phylogenetic tree from Bracken/Kraken2 reports by pruning GTDB/NCBI trees"
    )
    parser.add_argument(
        "--input_dir",
        "-i",
        required=True,
        help="Directory containing Bracken/Kraken2 report files",
    )
    parser.add_argument(
        "--bac_taxonomy",
        help="Path to GTDB bacterial taxonomy file",
    )
    parser.add_argument(
        "--ar_taxonomy",
        help="Path to GTDB archaeal taxonomy file",
    )
    parser.add_argument(
        "--bac_tree",
        help="Path to GTDB bacterial tree file",
    )
    parser.add_argument(
        "--ar_tree",
        help="Path to GTDB archaeal tree file",
    )
    parser.add_argument(
        "--out_prefix",
        "-o",
        default="output",
        help="Prefix for output files (default: output). Creates <prefix>.tree and <prefix>.otu.csv",
    )
    parser.add_argument(
        "--mode",
        "-m",
        choices=["bracken", "kraken2"],
        default="bracken",
        help="Input file format (default: bracken)",
    )
    parser.add_argument(
        "--keep-spaces",
        action="store_true",
        default=False,
        help="Keep spaces in species names (default: False)",
    )
    parser.add_argument(
        "--taxonomy",
        "-t",
        choices=["gtdb", "ncbi"],
        default="ncbi",
        help="Taxonomy source to use (default: ncbi).",
    )
    parser.add_argument(
        "--full-taxonomy",
        "-f",
        action="store_true",
        default=False,
        help="Include full taxonomy info in OTU table (default: False)",
    )
    args = parser.parse_args()
    if args.taxonomy == "gtdb":
        missing_args = []
        if not args.bac_taxonomy:
            missing_args.append("--bac_taxonomy")
        if not args.ar_taxonomy:
            missing_args.append("--ar_taxonomy")
        if not args.bac_tree:
            missing_args.append("--bac_tree")
        if not args.ar_tree:
            missing_args.append("--ar_tree")
        if missing_args:
            parser.error(
                "For GTDB taxonomy, the following arguments are required: "
                + ", ".join(missing_args)
            )
    return args


args = parse_args()


def main():
    # load data and create OTU table
    otu_table = pd.DataFrame()
    file_pattern = "*.breport" if args.mode == "bracken" else "*.k2report"

    matching_files = glob.glob(os.path.join(args.input_dir, file_pattern))
    if not matching_files:
        print(f"Error: No {args.mode} files found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    for f in matching_files:
        name = os.path.splitext(os.path.basename(f))[0]
        if args.mode == "bracken":
            sample_otu = read_bracken_report(f, taxonomy=args.taxonomy)
        else:
            sample_otu = read_kraken2_report(f, taxonomy=args.taxonomy)

        # Check if sample_otu is empty or doesn't have the required columns
        if sample_otu.empty or not all(
            col in sample_otu.columns
            for col in (
                ["taxid"] if args.taxonomy == "ncbi" else ["species", "abundance"]
            )
        ):
            print(
                f"Warning: Skipping {f} because it's empty or missing required columns.",
                file=sys.stderr,
            )
            continue

        sample_otu["sample"] = name
        otu_table = pd.concat([otu_table, sample_otu], ignore_index=True)

    if args.taxonomy == "ncbi":
        from ete3 import NCBITaxa

        ncbi = NCBITaxa()

        # Get unique taxids and translate to species names
        taxid_list = list(otu_table["taxid"].unique())
        translator = ncbi.get_taxid_translator(taxid_list)
        otu_table["species"] = otu_table["taxid"].map(
            lambda tid: translator.get(int(tid), str(tid))
        )

        if not args.keep_spaces:
            otu_table["species"] = otu_table["species"].str.replace(" ", "_")

        wide_otu_table = otu_table.pivot_table(
            index="species", columns="sample", values="abundance", aggfunc="sum"
        )
        if not args.full_taxonomy:
            wide_otu_table = wide_otu_table.reset_index()  # ensure species is a column
        wide_otu_table.to_csv(f"{args.out_prefix}.otu.csv", sep=",", index=False)

        if args.full_taxonomy:
            mapping = nc.build_lineage_mapping(ncbi, otu_table)
            wide_otu_table = wide_otu_table.reset_index()
            for col in ["domain", "phylum", "class", "order", "family", "genus"]:
                wide_otu_table[col] = wide_otu_table["species"].map(
                    lambda sp: mapping.get(sp, {}).get(col, "")
                )
            other_cols = [
                c
                for c in wide_otu_table.columns
                if c
                not in [
                    "species",
                    "domain",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                ]
            ]
            wide_otu_table = wide_otu_table[
                ["domain", "phylum", "class", "order", "family", "genus", "species"]
                + other_cols
            ]
            wide_otu_table.to_csv(f"{args.out_prefix}.otu.csv", sep=",", index=False)

        # Build tree using NCBI taxonomy and update leaf names
        tree = ncbi.get_topology(taxid_list)
        for leaf in tree.get_leaves():
            leaf.name = translator.get(int(leaf.name), leaf.name)
            if not args.keep_spaces:
                leaf.name = leaf.name.replace(" ", "_")
        tree.write(outfile=f"{args.out_prefix}.tree", format=1, quoted_node_names=False)
    else:  # GTDB
        wide_otu_table = otu_table.pivot_table(
            index="species", columns="sample", values="abundance", aggfunc="sum"
        )
        if args.full_taxonomy:
            wide_otu_table = wide_otu_table.reset_index()
            if not args.keep_spaces:
                wide_otu_table["species"] = wide_otu_table["species"].str.replace(
                    " ", "_"
                )
        else:
            wide_otu_table.index = wide_otu_table.index.str.replace(
                r"^s__", "", regex=True
            )
            if not args.keep_spaces:
                wide_otu_table.index = wide_otu_table.index.str.replace(" ", "_")
            wide_otu_table = (
                wide_otu_table.reset_index()
            )  # ensure species column is preserved

        if args.full_taxonomy:
            mapping = gt.build_taxonomy_mapping(
                [args.bac_taxonomy, args.ar_taxonomy], args.keep_spaces
            )
            for col in ["domain", "phylum", "class", "order", "family", "genus"]:
                wide_otu_table[col] = wide_otu_table["species"].map(
                    lambda sp: mapping.get(sp, {}).get(col, "")
                )
            other_cols = [
                c
                for c in wide_otu_table.columns
                if c
                not in [
                    "species",
                    "domain",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                ]
            ]
            wide_otu_table = wide_otu_table[
                ["domain", "phylum", "class", "order", "family", "genus", "species"]
                + other_cols
            ]

        wide_otu_table.to_csv(f"{args.out_prefix}.otu.csv", sep=",", index=False)

        # unique species (use original names to match taxonomy)
        species_set = set(otu_table["species"])

        # load taxonomies
        bac_species_to_genomes, bac_genome_to_species = gt.process_taxonomy(
            args.bac_taxonomy
        )
        ar_species_to_genomes, ar_genome_to_species = gt.process_taxonomy(
            args.ar_taxonomy
        )

        my_species = {"s__" + x if not x.startswith("s__") else x for x in species_set}

        bac_found, _ = gt.find_matching_species(my_species, bac_species_to_genomes)
        ar_found, _ = gt.find_matching_species(my_species, ar_species_to_genomes)

        bac_tree = gt.process_tree(args.bac_tree, bac_genome_to_species, bac_found)
        ar_tree = gt.process_tree(args.ar_tree, ar_genome_to_species, ar_found)

        combined_tree = Tree3()
        combined_tree.name = "root"
        combined_tree.add_child(bac_tree)
        combined_tree.add_child(ar_tree)
        combined_tree.write(
            outfile=f"{args.out_prefix}.tree", format=1, quoted_node_names=False
        )


if __name__ == "__main__":
    main()
