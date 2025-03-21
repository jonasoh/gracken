import pandas as pd


def read_bracken_report(file_path):
    """Extract species abundances from Bracken report."""
    df = pd.read_csv(file_path, sep="\t", header=None)
    filtered = df[(df[3] == "S") & (df[2] > 0)]
    result = pd.DataFrame()
    result["abundance"] = filtered[2]
    result["taxid"] = filtered[4] if df.shape[1] > 4 else ""
    result["species"] = filtered[5].str.strip() if df.shape[1] > 5 else ""
    return result


def read_kraken2_report(file_path):
    """Extract species abundances from Kraken2 report."""
    df = pd.read_csv(file_path, sep="\t", header=None)

    # Determine the format based on the number of columns
    num_cols = df.shape[1]

    # Different KRaken2 formats:
    # XXX: there might be further formats
    # %coverage, #reads, #reads direct, rank code, NCBI ID, name
    if num_cols == 6:
        rank_col = 3
        taxid_col = 4
        name_col = 5
    # %coverage, #reads, #reads direct, ?, ?, rank code, NCBI ID, name
    elif num_cols == 8:
        rank_col = 5
        taxid_col = 6
        name_col = 7
    else:
        raise ValueError(f"Unsupported Kraken2 report format: {num_cols} columns")

    filtered = df[(df[rank_col] == "S") & (df[1] > 0)]

    result = pd.DataFrame()
    result["abundance"] = filtered[1]
    result["taxid"] = filtered[taxid_col]
    result["species"] = filtered[name_col].str.strip()

    return result
