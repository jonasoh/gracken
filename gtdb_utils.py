# helper functions for GTDB-specific operations


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
        # Using lstrip to harmonize the leaf name for matching
        if leaf.name.lstrip("s__") in harmonized_keep:
            valid_species.append(leaf.name)
    if not valid_species:
        raise ValueError("No matching species found in tree")
    tree.prune(valid_species, preserve_branch_length=True)
    return tree
