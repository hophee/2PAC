# 2PAC

## Installation
Use `install.sh` to install most dependencies. You will also need to install the following R packages: `dplyr`, `readr`, `stringr`, `Biostrings`, `read.gb`.

## Usage
`./chopa.sh genome.fasta genome_annotation.tsv gene_name pTarget.db output_directory`

- genome_annotation.tsv accepted from [Bakta](https://github.com/oschwengers/bakta)
- gene_name may be in symbolic or `Locus Tag` format