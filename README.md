# Gracken

Gracken is a tool for creating phylogenetic trees from Bracken/Kraken2 reports. Both NCBI and GTDB taxonomies (see [Struo2](https://github.com/leylabmpi/Struo2)) are supported. It accomplishes this by pruning the NCBI/GTDB-provided tree using the species information from the report files. Gracken outputs a Newick-formatted tree and an OTU file with a format suitable for use with [phyloseq](https://joey711.github.io/phyloseq/).

Gracken was inspired by [this comment](https://github.com/jenniferlu717/KrakenTools/issues/46#issuecomment-2387744942) on the [KrakenTools](https://github.com/jenniferlu717/KrakenTools/) repository.

## Installation

If you have `uv` installed, you can run gracken directly using `uvx`:

```bash
uvx gracken --help
```

or install it using either `uv` or `pipx`:

```bash
# with uv
uv tool install gracken

# with pipx
pipx install gracken

# test that it works
gracken --help
```

## Usage

Gracken requires Bracken or Kraken2 report files. For GTDB taxonomy, you'll also need the taxonomy and tree files matching the database you used.

```
usage: gracken [-h] [--version] --input_dir INPUT_DIR [--bac_taxonomy BAC_TAXONOMY] [--ar_taxonomy AR_TAXONOMY]
               [--bac_tree BAC_TREE] [--ar_tree AR_TREE] [--out_prefix OUT_PREFIX] [--mode {bracken,kraken2}]
               [--keep-spaces] [--taxonomy {gtdb,ncbi}] [--full-taxonomy]

Creates a phylogenetic tree and OTU table from Bracken/Kraken2 reports by pruning GTDB/NCBI trees

options:
  -h, --help            show this help message and exit
  --version, -v         show program version and exit
  --input_dir INPUT_DIR, -i INPUT_DIR
                        directory containing Bracken/Kraken2 report files
  --bac_taxonomy BAC_TAXONOMY
                        path to GTDB bacterial taxonomy file
  --ar_taxonomy AR_TAXONOMY
                        path to GTDB archaeal taxonomy file
  --bac_tree BAC_TREE   path to GTDB bacterial tree file
  --ar_tree AR_TREE     path to GTDB archaeal tree file
  --out_prefix OUT_PREFIX, -o OUT_PREFIX
                        prefix for output files (default: output). Creates <prefix>.tree and <prefix>.otu.csv
  --mode {bracken,kraken2}, -m {bracken,kraken2}
                        input file format (default: bracken)
  --keep-spaces         keep spaces in species names (default: False)
  --taxonomy {gtdb,ncbi}, -t {gtdb,ncbi}
                        taxonomy source to use (default: ncbi).
  --full-taxonomy, -f   include full taxonomy info in OTU table (default: False)
```

### GTDB Taxonomy Example

First, download the taxonomy and tree files for the GTDB release corresponding to your Bracken reports from the [GTDB website](https://gtdb.ecogenomic.org/downloads). For example, if you used GTDB release 95, you would download these files: `bac120_taxonomy_r95.tsv`, `ar122_taxonomy_r95.tsv`, `bac120_r95.tree`, and `ar122_r95.tree`.

```bash
gracken --input_dir bracken_reports --bac_taxonomy bac120_taxonomy_r95.tsv --ar_taxonomy ar122_taxonomy_r95.tsv --bac_tree bac120_r95.tree --ar_tree ar122_r95.tree --out_prefix my_tree --taxonomy gtdb
```

This command will:

1.  Read Bracken reports from the `bracken_reports` directory.
2.  Use the specified GTDB taxonomy and tree files.
3.  Create two output files: `my_tree.tree` (the phylogenetic tree) and `my_tree.otu.csv` (the OTU table).

### NCBI Taxonomy Example

In NCBI mode, the taxdump is automatically downloaded. Simply run:

```bash
gracken --input_dir bracken_reports --out_prefix my_ncbi_tree --mode bracken
```

This command will:
1.  Read Bracken reports from the `bracken_reports` directory.
2.  Automatically download the necessary NCBI taxdump.
3.  Create two output files: `my_ncbi_tree.tree` (the phylogenetic tree) and `my_ncbi_tree.otu.csv` (the OTU table).

## Output Files

Gracken will output two files:

* `<prefix>.otu.csv`: A CSV file containing the OTU table. The first column contains the species names, and the remaining columns contain the counts for each species in each sample (smple names are derived from report file names). If you used the `--full-taxonomy` option, the first columns will contain the full taxonomy for each species.

* `<prefix>.tree`: A Newick-formatted tree file.

## Phyloseq Example

You can use the output OTU file and tree with [phyloseq](https://joey711.github.io/phyloseq/) and [ape](https://github.com/emmanuelparadis/ape). Here's an example R script which works regardless of whether you used `--full-taxonomy` or not:

```R
library(ape)
library(phyloseq)

# read the otu table and convert it to matrix, using the species column as row names
otu_tbl <- read.csv("output.otu.csv", stringsAsFactors = FALSE)
otu_mat <- as.matrix(otu_tbl[, (which(names(otu_tbl) == "species") + 1):ncol(otu_tbl)])
rownames(otu_mat) <- otu_tbl$species
otu_mat[is.na(otu_mat)] <- 0
otu_table_ps <- otu_table(otu_mat, taxa_are_rows = TRUE)

# load the tree
tree_ps <- phy_tree(read.tree('output.tree'))

# create the phyloseq object
ps <- phyloseq(otu_table_ps, tree_ps)
```
