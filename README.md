# Gracken

Gracken is a tool for creating phylogenetic trees from Bracken/Kraken2 reports. Both NCBI and GTDB taxonomies (see [Struo2](https://github.com/leylabmpi/Struo2)) are supported. It accomplishes this by pruning the NCBI/GTDB-provided tree using the species information from the report files. Gracken outputs a Newick-formatted tree and a OTU file with a format suitable for use with [phyloseq](https://joey711.github.io/phyloseq/).

Gracken was inspired by [this comment](https://github.com/jenniferlu717/KrakenTools/issues/46#issuecomment-2387744942) on the [KrakenTools](https://github.com/jenniferlu717/KrakenTools/) repository.

## Usage

Gracken requires Bracken or Kraken2 report files. For GTDB taxonomy, you'll also need the taxonomy and tree files matching the database you used.

```
usage: gracken.py  [-h] --input_dir INPUT_DIR [--bac_taxonomy BAC_TAXONOMY]
                   [--ar_taxonomy AR_TAXONOMY] [--bac_tree BAC_TREE]
                   [--ar_tree AR_TREE] [--out_prefix OUT_PREFIX]
                   [--mode {bracken,kraken2}] [--keep-spaces]
                   [--taxonomy {gtdb,ncbi}]
```

Options:
  -h, --help            Show this help message and exit.  
  --input_dir INPUT_DIR  Directory containing Bracken/Kraken2 report files.  
  --bac_taxonomy BAC_TAXONOMY  
                        Path to bacterial taxonomy file (e.g., bac120_taxonomy_r95.tsv).  
  --ar_taxonomy AR_TAXONOMY  
                        Path to archaeal taxonomy file (e.g., ar122_taxonomy_r95.tsv).  
  --bac_tree BAC_TREE   Path to bacterial tree file (e.g., bac120_r95.tree).  
  --ar_tree AR_TREE     Path to archaeal tree file (e.g., ar122_r95.tree).  
  --out_prefix OUT_PREFIX  
                        Prefix for output files (default: output). Creates <prefix>.tree and <prefix>.otu.csv.  
  --mode {bracken,kraken2}  
                        Input file format (default: bracken).  
  --keep-spaces  
                        Keep spaces in species names (default: False).  
  --taxonomy {gtdb,ncbi}  
                        Taxonomy source to use (default: ncbi).

### GTDB Taxonomy Example

First, download the taxonomy and tree files for the GTDB release corresponding to your Bracken reports from the [GTDB website](https://gtdb.ecogenomic.org/downloads). For example, if you used GTDB release 95, you would download these files: `bac120_taxonomy_r95.tsv`, `ar122_taxonomy_r95.tsv`, `bac120_r95.tree`, and `ar122_r95.tree`.

```bash
python gracken.py --input_dir bracken_reports --bac_taxonomy bac120_taxonomy_r95.tsv --ar_taxonomy ar122_taxonomy_r95.tsv --bac_tree bac120_r95.tree --ar_tree ar122_r95.tree --out_prefix my_tree --taxonomy gtdb
```

This command will:

1.  Read Bracken reports from the `bracken_reports` directory.
2.  Use the specified GTDB taxonomy and tree files.
3.  Create two output files: `my_tree.tree` (the phylogenetic tree) and `my_tree.otu.csv` (the OTU table).

### NCBI Taxonomy Example

In NCBI mode, the taxdump is automatically downloaded. Simply run:

```bash
python gracken.py --input_dir bracken_reports --out_prefix my_ncbi_tree --mode bracken
```

This command will:
1.  Read Bracken reports from the `bracken_reports` directory.
2.  Automatically download the necessary NCBI taxdump.
3.  Create two output files: `my_ncbi_tree.tree` (the phylogenetic tree) and `my_ncbi_tree.otu.csv` (the OTU table).

## Dependencies

*   Python 3
*   pandas
*   ete3

You can install the dependencies using pip:

```bash
pip install -r requirements.txt
```
