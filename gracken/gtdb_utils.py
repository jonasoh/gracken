from .vars import tax_cols


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
    """Checks the species_set against the full list of species."""
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
    from ete3 import Tree as Tree3

    tree = Tree3(tree_file, format=1, quoted_node_names=True)

    # replace tip names (genome ids) with species names from taxonomy
    for node in tree.traverse():
        if node.name in genome_to_species:
            node.name = genome_to_species[node.name]

    # Compute harmonized species to keep
    harmonized_keep = set(s.lstrip("s__") for s in species_to_keep)
    valid_species = []
    for leaf in tree.get_leaves():
        if leaf.name.lstrip("s__") in harmonized_keep:
            valid_species.append(leaf.name)

    if not valid_species:
        raise ValueError("No matching species found in tree")

    tree.prune(valid_species, preserve_branch_length=True)
    return tree


def build_taxonomy_mapping(tax_files, keep_spaces):
    """
    Build a mapping from GTDB taxonomy files.
    """
    gtdb_taxonomy = {}
    for tax_file in tax_files:
        with open(tax_file, "r") as tf:
            for line in tf:
                genome, tax_str = line.strip().split("\t")
                parts = tax_str.split(";")
                tax_dict = {}
                for part in parts:
                    for c in tax_cols:
                        if part.startswith(f"{c[0]}__"):
                            tax_dict[c] = part.replace(f"{c[0]}__", "")
                if "species" in tax_dict and not keep_spaces:
                    tax_dict["species"] = tax_dict["species"].replace(" ", "_")
                if "species" in tax_dict and tax_dict["species"] not in gtdb_taxonomy:
                    gtdb_taxonomy[tax_dict["species"]] = tax_dict
    return gtdb_taxonomy
