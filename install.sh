#!/usr/bin/env bash
set -uo pipefail

readonly PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly ENV_FILE="$PROJECT_DIR/env.yml"
readonly CONDA_ENV_NAME="oligo_design"
readonly CHOPCHOP_ENV_FILE="$PROJECT_DIR/env_chopchop.yml"
readonly CHOPCHOP_ENV_NAME="oligo_design_chopchop"
readonly VIENNARNA_ENV_FILE="$PROJECT_DIR/env_viennarna.yml"
readonly VIENNARNA_ENV_NAME="oligo_design_viennarna"
readonly CALL_PRIMER3_URL="https://gist.githubusercontent.com/IdoBar/5e78ae7a5cc7277a04b126ce6f595d6e/raw/45c60662f3479f41765bce839835c4988a7e5b36/callPrimer3.R"
readonly TEST_DIR="$PROJECT_DIR/test"
readonly MELTING_WRAPPER="$PROJECT_DIR/tools/melting-batch"
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
  local environment_name=$1
  local environment_file=$2
  [[ -f "$environment_file" ]] || {
    printf 'ERROR: Conda environment file not found: %s\n' \
      "$environment_file" >&2
    return 1
  }
  eval "$(conda shell.bash hook)" || return 1

  if conda env list | awk 'NR > 2 { print $1 }' | grep -Fqx "$environment_name"; then
    conda env update --name "$environment_name" --file "$environment_file" --prune
  else
    conda env create --file "$environment_file" --yes
  fi
}

install_melting_wrapper() {
  local environment_prefix
  environment_prefix="$(
    conda env list | awk -v name="$CONDA_ENV_NAME" '$1 == name { print $NF }'
  )"
  [[ -n "$environment_prefix" && -d "$environment_prefix/bin" ]] || {
    printf 'ERROR: cannot find Conda environment prefix for %s.\n' \
      "$CONDA_ENV_NAME" >&2
    return 1
  }
  [[ -x "$MELTING_WRAPPER" ]] || chmod +x "$MELTING_WRAPPER" || return 1
  ln -sfn "$MELTING_WRAPPER" "$environment_prefix/bin/melting-batch"
}

install_viennarna_commands() {
  local main_prefix
  local viennarna_prefix
  local command_name
  main_prefix="$(
    conda env list | awk -v name="$CONDA_ENV_NAME" '$1 == name { print $NF }'
  )"
  viennarna_prefix="$(
    conda env list | awk -v name="$VIENNARNA_ENV_NAME" \
      '$1 == name { print $NF }'
  )"
  [[ -n "$main_prefix" && -n "$viennarna_prefix" && \
    -d "$main_prefix/bin" ]] || {
    printf 'ERROR: cannot locate the main or ViennaRNA environment.\n' >&2
    return 1
  }
  for command_name in RNAfold RNAduplex RNAsubopt; do
    [[ -x "$viennarna_prefix/bin/$command_name" ]] || {
      printf 'ERROR: ViennaRNA command not found: %s\n' "$command_name" >&2
      return 1
    }
    ln -sfn "$viennarna_prefix/bin/$command_name" \
      "$main_prefix/bin/$command_name" || return 1
  done
}

install_chopchop_python() {
  local main_prefix
  local chopchop_prefix
  main_prefix="$(
    conda env list | awk -v name="$CONDA_ENV_NAME" '$1 == name { print $NF }'
  )"
  chopchop_prefix="$(
    conda env list | awk -v name="$CHOPCHOP_ENV_NAME" \
      '$1 == name { print $NF }'
  )"
  [[ -n "$main_prefix" && -n "$chopchop_prefix" && \
    -d "$main_prefix/bin" && -x "$chopchop_prefix/bin/python" ]] || {
    printf 'ERROR: cannot locate the isolated CHOPCHOP Python 2 runtime.\n' >&2
    return 1
  }
  ln -sfn "$chopchop_prefix/bin/python" "$main_prefix/bin/chopchop-python"
}

verify_primer_qc_dependencies() {
  conda run --no-capture-output --name "$CONDA_ENV_NAME" Rscript -e '
    oligoarray_data <- file.path(Sys.getenv("CONDA_PREFIX"), "share", "oligoarrayaux")
    if (file.exists(file.path(oligoarray_data, "stack.DAT"))) {
      Sys.setenv(UNAFOLDDAT = oligoarray_data)
    }
    stopifnot(
      requireNamespace("Biostrings", quietly = TRUE),
      requireNamespace("openPrimeR", quietly = TRUE),
      nzchar(Sys.which("hybrid-min")),
      nzchar(Sys.which("RNAfold")),
      nzchar(Sys.which("melting-batch")),
      nzchar(Sys.which("chopchop-python"))
    )
    settings <- openPrimeR::read_settings(system.file(
      "extdata", "settings", "C_Taq_PCR_high_stringency.xml",
      package = "openPrimeR"
    ))
    required <- c(
      "primer_length", "gc_ratio", "gc_clamp", "no_runs", "no_repeats",
      "melting_temp_range", "melting_temp_diff", "self_dimerization",
      "cross_dimerization", "secondary_structure", "primer_coverage"
    )
    stopifnot(all(required %in% names(openPrimeR::constraints(settings))))
  '
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
  run_step "Conda environment $CONDA_ENV_NAME" \
    install_conda_environment "$CONDA_ENV_NAME" "$ENV_FILE"
  run_step "Conda environment $CHOPCHOP_ENV_NAME" \
    install_conda_environment "$CHOPCHOP_ENV_NAME" "$CHOPCHOP_ENV_FILE"
  run_step "Conda environment $VIENNARNA_ENV_NAME" \
    install_conda_environment "$VIENNARNA_ENV_NAME" "$VIENNARNA_ENV_FILE"
  run_step "MELTING command wrapper" install_melting_wrapper
  run_step "isolated CHOPCHOP Python 2" install_chopchop_python
  run_step "isolated ViennaRNA commands" install_viennarna_commands
  run_step "Primer QC dependencies" verify_primer_qc_dependencies
else
  printf 'FAILED: Conda environment %s (conda was not found in PATH)\n' \
    "$CONDA_ENV_NAME" >&2
  FAILED_COMPONENTS+=("Conda environment $CONDA_ENV_NAME")
fi
run_step "primer3 source" clone_if_missing https://github.com/primer3-org/primer3.git primer3
run_step "primer3 build" make -C primer3/src
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
