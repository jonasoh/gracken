import pandas as pd


def read_bracken_report(file_path, taxonomy="ncbi"):
    """Extract species abundances from Bracken report."""
    df = pd.read_csv(file_path, sep="\t", header=None)
    if taxonomy == "ncbi":
        species_df = df[(df[3] == "S") & (df[2] > 0)]
        species_abundance = species_df[[4, 2]]
        species_abundance.columns = ["taxid", "abundance"]
    else:
        species_df = df[(df[3] == "S") & (df[2] > 0)]
        species_abundance = species_df[[5, 2]]
        species_abundance.columns = ["species", "abundance"]
        species_abundance["species"] = species_abundance["species"].str.strip()
    return species_abundance


def read_kraken2_report(file_path, taxonomy="ncbi"):
    """Extract species abundances from Kraken2 report."""
    df = pd.read_csv(file_path, sep="\t", header=None)
    if taxonomy == "ncbi":
        species_df = df[(df[5] == "S") & (df[1] > 0)]
        species_abundance = species_df[[6, 3]]
        species_abundance.columns = ["taxid", "abundance"]
    else:
        species_df = df[(df[5] == "S") & (df[1] > 0)]
        species_abundance = species_df[[7, 3]]
        species_abundance.columns = ["species", "abundance"]
        species_abundance["species"] = species_abundance["species"].str.strip()
    return species_abundance
