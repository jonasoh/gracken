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
    filtered = df[(df[5] == "S") & (df[1] > 0)]
    result = pd.DataFrame()
    result["abundance"] = filtered[3]
    result["taxid"] = filtered[6] if df.shape[1] > 6 else ""
    result["species"] = filtered[7].str.strip() if df.shape[1] > 7 else ""
    return result
