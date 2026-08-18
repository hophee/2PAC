#!/usr/bin/env bash
set -uo pipefail

readonly TEST_SCRIPT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")"
readonly CONDA_ENV_NAME="oligo_design"

if [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV_NAME" ]]; then
  command -v conda >/dev/null 2>&1 || {
    printf 'TEST FAILED: conda is not available in PATH\n' >&2
    exit 1
  }
  exec conda run --no-capture-output --name "$CONDA_ENV_NAME" \
    bash "$TEST_SCRIPT" "$@"
fi

readonly TEST_DIR="$(dirname -- "$TEST_SCRIPT")"
readonly PROJECT_DIR="$(cd -- "$TEST_DIR/.." && pwd)"
readonly OUTPUT_DIR="$TEST_DIR/test_output"
readonly WET_LAB_DIR="$OUTPUT_DIR/WetLab"
readonly TECH_REPORT_DIR="$OUTPUT_DIR/TechReport"
readonly SUMMARY="$TECH_REPORT_DIR/design_summary.tsv"
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

Rscript test/test_unit.R || fail "unit tests failed"

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
[[ ! -e "$OUTPUT_DIR/design_summary.tsv" ]] ||
  fail "design_summary.tsv was written outside TechReport"
[[ -s "$TECH_REPORT_DIR/run_parameters.tsv" ]] ||
  fail "run_parameters.tsv was not created"
[[ -s "$TECH_REPORT_DIR/chopchop_config.json" ]] ||
  fail "CHOPCHOP configuration snapshot was not created"
[[ -s "$TECH_REPORT_DIR/genome_indexes/MG1655.2bit" ]] ||
  fail "2bit genome index was not written to TechReport"
[[ -s "$TECH_REPORT_DIR/genome_indexes/MG1655.1.ebwt" ]] ||
  fail "Bowtie genome index was not written to TechReport"
grep -q $'^genome_file\t' "$TECH_REPORT_DIR/run_parameters.tsv" ||
  fail "run_parameters.tsv lacks the genome input path"
grep -q $'^n20_arm_min_distance_nt\t40$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
  fail "run_parameters.tsv lacks the N20-to-arm distance"
grep -q $'^primer3_buffer_divalent_salt_mm\t1.5$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
  fail "run_parameters.tsv lacks Primer3 buffer data"

for gene in recA pta hupB; do
  awk -F '\t' -v gene="$gene" 'NR > 1 && $1 == gene && $4 == "ok" { found = 1 } END { exit !found }' "$SUMMARY" ||
    fail "$gene is not marked as successful in design_summary.tsv"
  awk -F '\t' -v gene="$gene" 'NR > 1 && $1 == gene && $7 != "" && $7 != "NA" { found = 1 } END { exit !found }' "$SUMMARY" ||
    fail "$gene has no WetLab path in design_summary.tsv"

  target_dir="$TECH_REPORT_DIR/${gene,,}_results"
  wet_target_dir="$WET_LAB_DIR/${gene,,}_results"
  for result in all_primers.fasta edited_genome.fasta report.tsv n20_table.tsv design.log; do
    [[ -s "$target_dir/$result" ]] || fail "missing or empty result for $gene: $result"
  done
  for result in final_sequences.fasta final_sequences.txt wet_lab_report.txt; do
    [[ -s "$wet_target_dir/$result" ]] ||
      fail "missing or empty WetLab result for $gene: $result"
  done
  [[ ! -e "$wet_target_dir/design.log" ]] ||
    fail "WetLab unexpectedly contains a technical design.log for $gene"
  [[ ! -e "$wet_target_dir/n20_table.tsv" ]] ||
    fail "WetLab unexpectedly contains a CHOPCHOP table for $gene"
  fasta_records="$(grep -c '^>' "$wet_target_dir/final_sequences.fasta")"
  txt_records="$(awk 'NR > 1 { count++ } END { print count + 0 }' "$wet_target_dir/final_sequences.txt")"
  [[ "$fasta_records" -eq "$txt_records" ]] ||
    fail "WetLab FASTA and TXT sequence sets differ for $gene"
  grep -q 'Неуспешная вставка' "$wet_target_dir/wet_lab_report.txt" ||
    fail "WetLab report lacks the unsuccessful-insertion PCR size for $gene"
  grep -q 'Успешная вставка' "$wet_target_dir/wet_lab_report.txt" ||
    fail "WetLab report lacks the successful-insertion PCR size for $gene"
done

printf 'TEST PASSED: the MG1655 recA, pta and hupB test runs completed and expected output files are valid.\n'
