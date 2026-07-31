#!/usr/bin/env bash
set -uo pipefail

readonly PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly ENV_FILE="$PROJECT_DIR/env.yml"
readonly CONDA_ENV_NAME="oligo_design"
readonly CALL_PRIMER3_URL="https://gist.githubusercontent.com/IdoBar/5e78ae7a5cc7277a04b126ce6f595d6e/raw/45c60662f3479f41765bce839835c4988a7e5b36/callPrimer3.R"
readonly TEST_DIR="$PROJECT_DIR/test"
readonly MG1655_REFSEQ_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2"
readonly MG1655_GENOME_URL="$MG1655_REFSEQ_URL/GCF_000005845.2_ASM584v2_genomic.fna.gz"
readonly MG1655_GFF_URL="$MG1655_REFSEQ_URL/GCF_000005845.2_ASM584v2_genomic.gff.gz"

declare -a FAILED_COMPONENTS=()

run_step() {
  local component=$1
  shift

  if "$@"; then
    printf 'OK: %s\n' "$component"
  else
    printf 'FAILED: %s\n' "$component" >&2
    FAILED_COMPONENTS+=("$component")
  fi
}

clone_if_missing() {
  local repository=$1
  local destination=$2

  if [[ -e "$destination" && ! -d "$destination/.git" ]]; then
    printf 'ERROR: %s already exists but is not a Git checkout.\n' "$destination" >&2
    return 1
  fi
  if [[ ! -e "$destination" ]]; then
    GIT_TERMINAL_PROMPT=0 git clone --depth 1 "$repository" "$destination" || return 1
  fi
}

install_conda_environment() {
  [[ -f "$ENV_FILE" ]] || {
    printf 'ERROR: Conda environment file not found: %s\n' "$ENV_FILE" >&2
    return 1
  }
  eval "$(conda shell.bash hook)" || return 1

  if conda env list | awk 'NR > 2 { print $1 }' | grep -Fqx "$CONDA_ENV_NAME"; then
    conda env update --name "$CONDA_ENV_NAME" --file "$ENV_FILE" --prune
  else
    conda env create --file "$ENV_FILE" --yes
  fi
}

download_call_primer3() {
  local temporary_file
  temporary_file="$(mktemp "$PROJECT_DIR/.callPrimer3.R.XXXXXX")" || return 1

  if ! curl --fail --location --silent --show-error "$CALL_PRIMER3_URL" --output "$temporary_file"; then
    rm -f "$temporary_file"
    return 1
  fi
  mv "$temporary_file" "$PROJECT_DIR/callPrimer3.R"
}

download_gzip_if_missing() {
  local url=$1
  local destination=$2
  local temporary_file

  [[ -s "$destination" ]] && return 0
  temporary_file="$(mktemp "$TEST_DIR/.download.XXXXXX.gz")" || return 1
  if ! curl --fail --location --silent --show-error "$url" --output "$temporary_file"; then
    rm -f "$temporary_file"
    return 1
  fi
  if ! gzip -t "$temporary_file" || ! gzip -dc "$temporary_file" > "$destination"; then
    rm -f "$temporary_file" "$destination"
    return 1
  fi
  rm -f "$temporary_file"
}

prepare_test_data() {
  mkdir -p "$TEST_DIR" || return 1
  download_gzip_if_missing "$MG1655_GENOME_URL" "$TEST_DIR/MG1655.fna" || return 1
  download_gzip_if_missing "$MG1655_GFF_URL" "$TEST_DIR/MG1655.gff" || return 1

  if [[ ! -f "$TEST_DIR/test_run.sh" ]]; then
    printf 'ERROR: test runner was not found: %s\n' "$TEST_DIR/test_run.sh" >&2
    return 1
  fi
  chmod +x "$TEST_DIR/test_run.sh"
}

patch_chopchop_python2_pandas_compatibility() {
  local chopchop_file="$PROJECT_DIR/chopchop/chopchop.py"
  local old_call='args.backbone, args.replace5P, args.maxOffTargets, countMM, args.PAM,'
  local new_call='args.backbone, args.replace5P, int(args.maxOffTargets), countMM, args.PAM,'

  [[ -f "$chopchop_file" ]] || return 1
  if grep -Fq "$new_call" "$chopchop_file"; then
    return 0
  fi
  if ! grep -Fq "$old_call" "$chopchop_file"; then
    printf 'ERROR: CHOPCHOP source line for maxOffTargets patch was not found.\n' >&2
    return 1
  fi
  sed -i "s/${old_call}/${new_call}/" "$chopchop_file"
}

cd "$PROJECT_DIR"
if command -v conda >/dev/null; then
  run_step "Conda environment $CONDA_ENV_NAME" install_conda_environment
else
  printf 'FAILED: Conda environment %s (conda was not found in PATH)\n' \
    "$CONDA_ENV_NAME" >&2
  FAILED_COMPONENTS+=("Conda environment $CONDA_ENV_NAME")
fi
run_step "primer3 source" clone_if_missing https://github.com/primer3-org/primer3.git primer3
run_step "primer3 build" make -C primer3/src
run_step "virtualPCR source" clone_if_missing https://github.com/rkalendar/virtualPCR.git virtualPCR
run_step "CHOPCHOP source" clone_if_missing https://github.com/JokingHero/chopchop.git chopchop
run_step "CHOPCHOP Python 2/pandas compatibility" patch_chopchop_python2_pandas_compatibility
run_step "callPrimer3.R" download_call_primer3
run_step "MG1655 test data and test runner" prepare_test_data

if (( ${#FAILED_COMPONENTS[@]} == 0 )); then
  printf 'Installation completed successfully.\n'
else
  printf '\nInstallation completed with missing components:\n' >&2
  printf '  - %s\n' "${FAILED_COMPONENTS[@]}" >&2
  exit 1
fi
