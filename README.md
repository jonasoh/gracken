# Gracken

Gracken is a simple tool for creating a phylogenetic tree from Bracken/Kraken2 reports when using GTDB taxonomy (see [Struo2](https://github.com/leylabmpi/Struo2)). It works by pruning the GTDB tree using the species information from the report files. This makes it possible to do e.g., [phyloseq](https://joey711.github.io/phyloseq/) analyses directly on these reports.

Gracken was inspired by [this comment](https://github.com/jenniferlu717/KrakenTools/issues/46#issuecomment-2387744942) on the [KrakenTools](https://github.com/jenniferlu717/KrakenTools/) repository.

## Usage

Gracken requires Bracken or Kraken2 report files, as well as the matching GTDB taxonomy and tree files.

```
usage: gracken.py [-h] --input_dir INPUT_DIR --bac_taxonomy BAC_TAXONOMY
                   --ar_taxonomy AR_TAXONOMY --bac_tree BAC_TREE --ar_tree
                   AR_TREE [--out_prefix OUT_PREFIX] [--mode {bracken,kraken2}]
                   [--dont-replace-spaces]

Create a phylogenetic tree using Bracken/Kraken2 reports by pruning the GTDB tree

options:
  -h, --help            show this help message and exit
  --input_dir INPUT_DIR
                        Directory containing Bracken report files
  --bac_taxonomy BAC_TAXONOMY
                        Path to bacterial taxonomy file (e.g.,
                        bac120_taxonomy_r95.tsv)
  --ar_taxonomy AR_TAXONOMY
                        Path to archaeal taxonomy file (e.g.,
                        ar122_taxonomy_r95.tsv)
  --bac_tree BAC_TREE   Path to bacterial tree file (e.g., bac120_r95.tree)
  --ar_tree AR_TREE     Path to archaeal tree file (e.g., ar122_r95.tree)
  --out_prefix OUT_PREFIX
                        Prefix for output files (default: output). Creates
                        <prefix>.tree and <prefix>.otu.csv
  --mode {bracken,kraken2}
                        Input file format (default: bracken)
  --dont-replace-spaces
                        Keep spaces in species names (default: False)
```

### Example

First, download the taxonomy and tree files matching the GTDB release used for your Bracken reports from the [GTDB website](https://gtdb.ecogenomic.org/downloads).

```bash
python gracken.py --input_dir bracken_reports --bac_taxonomy bac120_taxonomy_r95.tsv --ar_taxonomy ar122_taxonomy_r95.tsv --bac_tree bac120_r95.tree --ar_tree ar122_r95.tree --out_prefix my_tree
```

This command will:

1.  Read Bracken reports from the `bracken_reports` directory.
2.  Use the specified GTDB taxonomy and tree files.
3.  Create two output files: `my_tree.tree` (the phylogenetic tree) and `my_tree.otu.csv` (the OTU table).

## Dependencies

*   Python 3
*   pandas
*   ete3

You can install the dependencies using pip:

```bash
pip install -r requirements.txt
```
