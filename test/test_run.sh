#!/usr/bin/env bash
set -uo pipefail

readonly TEST_SCRIPT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")"

if [[ "${CONDA_DEFAULT_ENV:-}" != "chopchop" ]]; then
  command -v conda >/dev/null 2>&1 || {
    printf 'TEST FAILED: conda is not available in PATH\n' >&2
    exit 1
  }
  exec conda run --no-capture-output --name chopchop bash "$TEST_SCRIPT" "$@"
fi

readonly TEST_DIR="$(dirname -- "$TEST_SCRIPT")"
readonly PROJECT_DIR="$(cd -- "$TEST_DIR/.." && pwd)"
readonly OUTPUT_DIR="$TEST_DIR/test_output"
readonly SUMMARY="$OUTPUT_DIR/design_summary.tsv"
readonly PLASMID="$TEST_DIR/test_target_plasmid.fasta"
readonly TEST_GENES="recA,pta,hupB"

fail() {
  printf 'TEST FAILED: %s\n' "$1" >&2
  exit 1
}

[[ -s "$TEST_DIR/MG1655.fna" ]] || fail "MG1655.fna is missing or empty"
[[ -s "$TEST_DIR/MG1655.gff" ]] || fail "MG1655.gff is missing or empty"

cat > "$PLASMID" <<'EOF'
>synthetic_test_target_plasmid
ACGTACGATCGATGCTAGCTACGATCGTACGATCGATGCTAGCTACGATCGTACGATCGATGCTAGCTACGATCGT
EOF

rm -rf "$OUTPUT_DIR"
cd "$PROJECT_DIR" || fail "cannot enter project directory"

if ! Rscript oligo_designer.R \
  --genome "$TEST_DIR/MG1655.fna" \
  --genome-annotation "$TEST_DIR/MG1655.gff" \
  --annotation-format gff \
  --target-plasmid "$PLASMID" \
  --output-dir "$OUTPUT_DIR" \
  --cds "$TEST_GENES"; then
  fail "oligo_designer.R returned a non-zero exit code"
fi

[[ -s "$SUMMARY" ]] || fail "design_summary.tsv was not created"

for gene in recA pta hupB; do
  awk -F '\t' -v gene="$gene" 'NR > 1 && $1 == gene && $4 == "ok" { found = 1 } END { exit !found }' "$SUMMARY" ||
    fail "$gene is not marked as successful in design_summary.tsv"

  target_dir="$OUTPUT_DIR/${gene,,}_results"
  for result in all_primers.fasta edited_genome.fasta report.tsv n20_table.tsv; do
    [[ -s "$target_dir/$result" ]] || fail "missing or empty result for $gene: $result"
  done
done

printf 'TEST PASSED: the MG1655 recA, pta and hupB test runs completed and expected output files are valid.\n'
