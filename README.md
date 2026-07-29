# 2PAC

## Installation
Use `./install.sh` to install most dependencies. The script is safe to rerun: it updates the `chopchop` Conda environment and reuses existing Git checkouts. It stops at the first failed command.

The CHOPCHOP upstream repository is cloned from [GitHub](https://github.com/JokingHero/chopchop).

You will also need to install the following R packages: `argparser`, `dplyr`, `readr`, `Biostrings`, `ape`.

## Usage
```bash
Rscript oligo_designer.R \
  --genome genome.fasta --genome-annotation genome_annotation.tsv \
  --target-plasmid pTarget.fasta --output-dir output_directory \
  --cds gene1,gene2 --ncrna rna1,rna2
```

- `genome_annotation.tsv` accepted from [Bakta](https://github.com/oschwengers/bakta)
  (`--annotation-format bakta`) or GFF/GFF3 (`--annotation-format gff`, the default). Target
  is first looked up in `locus_tag`, then in `gene`.
- Gene names may be symbolic names or `Locus Tag` values. `--cds` and `--ncrna`
  accept comma-separated or space-separated lists.

## Output
- `all_primers.fasta`: set of designed primers
- `edited_genome.fasta` and `pTargets.fasta`: genome and 3 pTargets sequences after editing
- `n20_table.tsv`: chopchop.py output
- `offtarget_check.txt`: virtualPCR report
- `report.txt`: report on melting temperature and PCR product lengths
- `PCR_product_sequnces.fasta`: sequence of all PCR products
