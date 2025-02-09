from .vars import tax_cols


def build_lineage_mapping(ncbi, otu_table):
    """
    Build a mapping from species to taxonomy info using NCBI lineage.
    """
    species_taxid = otu_table.groupby("species")["taxid"].first().to_dict()
    lineage_mapping = {}
    for sp, tid in species_taxid.items():
        lineage = ncbi.get_lineage(int(tid))
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        taxonomy = {c: "" for c in tax_cols}
        for rank in taxonomy:
            if rank == "species":
                continue
            for t in lineage:
                actual_rank = ranks.get(t)
                if rank == "domain" and actual_rank == "superkingdom":
                    taxonomy[rank] = names.get(t, "")
                    break
                elif actual_rank == rank:
                    taxonomy[rank] = names.get(t, "")
                    break
        lineage_mapping[sp] = taxonomy
    return lineage_mapping
