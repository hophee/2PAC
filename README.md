# 2PAC

## Installation
Use `./install.sh` to install most dependencies. The script is safe to rerun: it updates the `chopchop` Conda environment and reuses existing Git checkouts. It stops at the first failed command.

The CHOPCHOP upstream repository is cloned from Bitbucket. If your account cannot access it, obtain a compatible checkout and either place it in `chopchop/` or set `CHOPCHOP_REPOSITORY` to an accessible Git URL before running the installer.

You will also need to install the following R packages: `dplyr`, `readr`, `stringr`, `Biostrings`, `read.gb`.

## Usage
`./chopa.sh genome.fasta genome_annotation.tsv gene_name pTarget.gb output_directory`

- genome_annotation.tsv accepted from [Bakta](https://github.com/oschwengers/bakta)
- gene_name may be in symbolic or `Locus Tag` format

## Output
- `all_primers.fasta`: set of designed primers
- `edited_genome.fasta` and `pTargets.fasta`: genome and 3 pTargets sequences after editing
- `n20_table.tsv`: chopchop.py output
- `offtarget_check.txt`: virtualPCR report
- `report.txt`: report on melting temperature and PCR product lengths
- `PCR_product_sequnces.fasta`: sequence of all PCR products
