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

The output root is split by audience:

```text
results/
├── WetLab/
│   └── <target>_results/
│       ├── final_sequences.fasta
│       ├── final_sequences.txt
│       └── wet_lab_report.txt
└── TechReport/
    ├── design_summary.tsv
    ├── run_parameters.tsv
    └── <target>_results/
        └── ... pipeline outputs, tool reports, and logs
```

`final_sequences.fasta` and `final_sequences.txt` contain the same complete
oligonucleotide set: forward N20 oligos, the common reverse sgRNA oligo, four
homology-arm primers, and two screening primers. `wet_lab_report.txt` lists
the annealing-region Tm values and the expected screening PCR product sizes
for the unedited (unsuccessful insertion) and edited (successful insertion)
alleles. A target-specific WetLab directory is written only after all design
and virtualPCR stages succeed.

`TechReport/` contains everything needed for traceability: CHOPCHOP, Primer3,
and virtualPCR products, the edited-genome model, intermediate FASTA/TSV
files, genome indexes, a CHOPCHOP configuration snapshot, `report.tsv`, and
`design.log`. `run_parameters.tsv` records input file paths, target lists, N20
and homology-arm settings, external tool paths, and the Primer3 buffer
concentrations used for Tm calculation (50 mM monovalent salt, 1.5 mM
divalent salt, 0.6 mM dNTP, and 50 nM DNA).

If a target fails, its `TechReport/<target>_results/` directory also contains
`error.txt` with the gene, design class, failed stage, and reason. Other
targets continue to run and the failure is recorded in
`TechReport/design_summary.tsv`; no new WetLab result is created for that
failed target.

The implementation still supports one complete linear genome contig only.
