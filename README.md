# 2PAC

2PAC designs CRISPR-Cas9 N20 sequences, homology-arm primers, screening
primers, and an edited-genome model for bacterial CDS and ncRNA targets.

## Installation

Run `./install.sh` to install the external tools and create the `oligo_design`
Conda environment. 
The required R packages are `argparser`, `dplyr`, `readr`, `Biostrings`, `ape`, and `janitor`.

## Usage

```bash
Rscript oligo_designer.R \
  --genome genome.fasta \
  --genome-annotation genome.gff \
  --annotation-format gff \
  --target-plasmid pTarget.fasta \
  --output-dir results \
  --cds gene1,gene2 \
  --ncrna rna1,rna2
```

Targets are first matched by `locus_tag`, then by `gene`. `--cds` and
`--ncrna` accept comma-separated or space-separated values. Bakta TSV and
GFF/GFF3 annotations are supported.

### Design parameters

| Argument | Default | Meaning |
|---|---:|---|
| `--n20-mn` | `1` | Number of N20 sequences required in the final set |
| `--n20-strands` | `random` | For sets larger than one: `plus`, `minus`, `both`, or `random` |
| `--n20-offtarget` | `0` | Maximum values for `MM0,MM1,...`, separated by commas |
| `--cds-fs` | `FALSE` | Require the deleted CDS length to be divisible by three |
| `--ncrna-fs` | `FALSE` | Require the deleted ncRNA length to be divisible by three |
| `--left-arm-min/opt/max` | `300/350/400` | Hard minimum, preferred initial, and hard maximum left-arm lengths |
| `--right-arm-min/opt/max` | `400/450/500` | Hard minimum, preferred initial, and hard maximum right-arm lengths |
| `--n20-arm-min-distance` | `40` | Minimum number of bases between every selected N20 and each homology arm |

`random` means that either strand composition is accepted; candidates remain
ranked deterministically. The number of `--n20-offtarget` values must not
exceed the number of `MM*` columns emitted by CHOPCHOP. A candidate is retained
only when every supplied threshold is met.

## Output

The output root contains `design_summary.tsv`. Each target has a
`<target>_results/` directory containing, among other intermediate files:

- `selected_n20_table.tsv` and `selected_n20.fasta`;
- `homology_arms.fasta`;
- `all_primers.fasta`;
- `primers_without_service_sequences.fasta`;
- `screening_primers.fasta`;
- `edited_genome.fasta`;
- `report.tsv`;
- virtualPCR reports;
- `design.log`.

If a target fails, its directory also contains `error.txt` with the gene,
design class, failed stage, and reason. Other targets continue to run and the
failure is recorded in `design_summary.tsv`.

The implementation still supports one complete linear genome contig only.
See `TECHNICAL_SPECIFICATION.md` for the full current contract and known
limitations.
