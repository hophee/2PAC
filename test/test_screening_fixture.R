#!/usr/bin/env Rscript

source("oligo_designer.R")

assert_true <- function(value, message) {
  if (!isTRUE(value)) {
    stop(message, call. = FALSE)
  }
}

local({
  target_dir <- tempfile("2pac-screening-fixture-")
  dir.create(target_dir)
  on.exit(unlink(target_dir, recursive = TRUE), add = TRUE)

  original_call_primer3 <- callPrimer3
  original_evaluate_candidate <- evaluate_candidate_reaction
  on.exit(
    assign("callPrimer3", original_call_primer3, envir = .GlobalEnv),
    add = TRUE
  )
  on.exit(
    assign(
      "evaluate_candidate_reaction",
      original_evaluate_candidate,
      envir = .GlobalEnv
    ),
    add = TRUE
  )

  assign(
    "callPrimer3",
    function(...) {
      data.frame(
        PRIMER_LEFT_SEQUENCE = c("ACGTCGATCGTAGCTACGTA", "GCTAGTCGATGCTACGTAGC"),
        PRIMER_RIGHT_SEQUENCE = c("TGCATCGATGCTAGTCGTAC", "CGTACGATCGTAGCATCGAC"),
        PRIMER_LEFT_pos = c(100L, 110L),
        PRIMER_RIGHT_pos = c(1099L, 1109L),
        PRIMER_LEFT_TM = c(62.5, 62.7),
        PRIMER_RIGHT_TM = c(62.6, 62.8),
        PRIMER_PAIR_PRODUCT_SIZE = c(1000L, 1000L),
        PRIMER_PAIR_PENALTY = c(0.1, 0.2),
        stringsAsFactors = FALSE
      )
    },
    envir = .GlobalEnv
  )

  assign(
    "evaluate_candidate_reaction",
    function(input, primer_row, full_forward, full_reverse, reaction, pair_id,
             trace) {
      passed <- identical(pair_id, "screening_02")
      specificity <- list(
        passed = passed,
        binding_sites = data.frame(
          primer_id = c("forward", "reverse"),
          reference_id = input$genome_reference_id,
          reference_type = "genome",
          start = c(primer_row$genome_start[[1]], primer_row$genome_end[[1]] - 19L),
          end = c(primer_row$genome_start[[1]] + 19L, primer_row$genome_end[[1]]),
          strand = c("+", "-"),
          mismatches = 0L,
          mismatch_positions = "",
          mismatches_3p = 0L,
          stringsAsFactors = FALSE
        ),
        amplicons = data.frame(
          reference_id = input$genome_reference_id,
          reference_type = "genome",
          start = primer_row$genome_start[[1]],
          end = primer_row$genome_end[[1]],
          product_size = 1000L,
          sequence = paste(rep("A", 1000L), collapse = ""),
          intended = passed,
          off_target = !passed,
          invalid_size = FALSE,
          circular_wrap = FALSE,
          rejection_reason = if (passed) "" else "off_target",
          stringsAsFactors = FALSE
        ),
        n_high_risk_offtarget_products = as.integer(!passed),
        n_all_offtarget_products = as.integer(!passed),
        n_perfect_3p_offtarget_sites = as.integer(!passed),
        rejection_reason = if (passed) "" else "off_target_products=1"
      )
      openprimer <- if (passed) {
        list(
          passed = TRUE,
          metrics = data.frame(
            reaction = reaction,
            constraints_passed = TRUE,
            penalty = 1,
            stringsAsFactors = FALSE
          ),
          failed_soft_constraints = 0L,
          penalty = 1,
          max_dimer_risk = 0,
          abs_tm_diff = 0.1,
          rejection_reason = ""
        )
      } else NULL
      append_primer_qc_trace(
        trace,
        reaction,
        pair_id,
        specificity = specificity,
        openprimer = openprimer
      )
      list(
        passed = passed,
        specificity = specificity,
        openprimer = openprimer,
        rejection_reason = specificity$rejection_reason
      )
    },
    envir = .GlobalEnv
  )

  trace <- new_primer_qc_trace()
  trace$ranking[[1]] <- data.frame(
    pair_id = "homology_fixture",
    reaction = "homology_combination",
    primer3_index = 1L,
    structure_passed = TRUE,
    specificity_passed = TRUE,
    openprimer_passed = TRUE,
    n_high_risk_offtarget_products = 0L,
    n_all_offtarget_products = 0L,
    n_perfect_3p_offtarget_sites = 0L,
    openprimer_failed_soft_constraints = 0L,
    openprimer_penalty = 1,
    max_dimer_risk = 0,
    abs_tm_diff = 0.1,
    primer3_pair_penalty = 0.1,
    deleted_nt = 200L,
    selected = TRUE,
    rejection_reason = "",
    stringsAsFactors = FALSE
  )

  pair <- data.frame(
    PRIMER_LEFT_SEQUENCE = c("ATGCGTACGATCGTACGTAG", "GATCGTAGCTAGTCGATGCA"),
    PRIMER_RIGHT_SEQUENCE = c("CTACGTACGATCGTACGCAT", "TGCATCGACTAGCTACGATC"),
    PRIMER_LEFT_pos = c(1L, 1L),
    PRIMER_RIGHT_pos = c(300L, 400L),
    PRIMER_LEFT_len = c(20L, 20L),
    PRIMER_RIGHT_len = c(20L, 20L),
    PRIMER_LEFT_TM = c(60, 60),
    PRIMER_RIGHT_TM = c(60, 60),
    genome_start = c(201L, 701L),
    genome_end = c(500L, 1100L),
    stringsAsFactors = FALSE
  )
  input <- list(
    genome = DNAString(paste(rep("ACGT", 500L), collapse = "")),
    genome_reference_id = "fixture_genome",
    parameters = list(
      cds_fs = FALSE,
      ncrna_fs = FALSE,
      n20_offtarget = 0L,
      n20_arm_min_distance = 40L
    ),
    tools = list(primer3 = "unused", primer3_config = "unused")
  )
  arms <- list(
    pair = pair,
    ticks = c(201L, 500L, 701L, 1100L),
    left = input$genome[201L:700L],
    right = input$genome[701L:1200L],
    selected_pair_id = "homology_fixture",
    primer_qc_trace = trace
  )
  selected <- list(table = data.frame(
    target_sequence = "ACGTACGTACGTACGTACGTAGG",
    strand = "+",
    n20_start = 560L,
    n20_end = 579L,
    stringsAsFactors = FALSE
  ))
  feature <- list(
    display_name = "fixture",
    query_name = "fixture",
    strand = "+"
  )
  log_path <- file.path(target_dir, "design.log")

  result <- write_design_outputs(
    input,
    feature,
    selected,
    arms,
    "cds",
    target_dir,
    log_path
  )
  write_primer_qc_trace(result$primer_qc_trace, target_dir)
  wet_lab_dir <- file.path(target_dir, "WetLab")
  write_wet_lab_outputs(
    wet_lab_dir,
    feature,
    "cds",
    result$wet_lab$sequences,
    result$wet_lab$sequence_purposes,
    result$wet_lab$primer_metrics,
    result$wet_lab$screening_product_sizes,
    result$wet_lab$n20_distances,
    result$wet_lab$screening_qc
  )

  assert_true(
    identical(result$screening_pair_id, "screening_02"),
    "The second screening pair was not selected"
  )
  screening <- as.character(readDNAStringSet(result$screening_path))
  assert_true(
    identical(
      unname(screening),
      c("GCTAGTCGATGCTACGTAGC", "CGTACGATCGTAGCATCGAC")
    ),
    "Screening output does not contain the selected Primer3 row"
  )
  ranking <- read.delim(
    file.path(target_dir, "primer_pair_ranking.tsv"),
    check.names = FALSE
  )
  screening_ranking <- ranking[ranking$reaction == "scrF_scrR", , drop = FALSE]
  assert_true(
    nrow(screening_ranking) == 2L &&
      !screening_ranking$selected[[1]] &&
      screening_ranking$selected[[2]],
    "Screening ranking did not reject row 1 and select row 2"
  )
  amplicons <- read.delim(
    file.path(target_dir, "primer_amplicons.tsv"),
    check.names = FALSE
  )
  assert_true(
    any(
      amplicons$pair_id == "screening_02" &
        amplicons$reaction == "scrF_scrR" &
        amplicons$intended &
        !amplicons$off_target
    ),
    "Selected screening pair lacks its intended amplicon"
  )
  log_lines <- readLines(log_path)
  assert_true(
    any(grepl("primer_qc\\tTRY\\tpair_id=screening_01", log_lines)) &&
      any(grepl("primer_qc\\tREJECTED\\tpair_id=screening_01", log_lines)) &&
      any(grepl("primer_qc\\tOK\\tpair_id=screening_02", log_lines)),
    "Screening TRY/REJECTED/OK events are incomplete"
  )
  report <- readLines(file.path(target_dir, "report.tsv"))
  assert_true(
    any(report == "screening_pair_id\tscreening_02"),
    "report.tsv lacks the selected screening pair ID"
  )
  assert_true(
    result$wet_lab$n20_distances$left_arm_distance_bp[[1]] == 59L &&
      result$wet_lab$n20_distances$right_arm_distance_bp[[1]] == 121L,
    "WetLab data lacks the selected N20-to-arm distances"
  )
  assert_true(
    result$wet_lab$screening_qc$offtarget_products == 0L &&
      any(result$wet_lab$screening_qc$openprimer_metrics$metric ==
        "Все обязательные ограничения пройдены"),
    "WetLab data lacks screening QC"
  )
  wet_lab_report <- readLines(
    file.path(wet_lab_dir, "wet_lab_report.txt"),
    encoding = "UTF-8"
  )
  assert_true(
    any(grepl("N20_1.*59.*121", wet_lab_report)) &&
      any(grepl("Оффтаргетные ПЦР-продукты, всего.*0", wet_lab_report)) &&
      any(grepl("Все обязательные ограничения пройдены.*пройдено", wet_lab_report)),
    "Selected screening QC was not written to the WetLab report"
  )
})

message("Screening integration fixture passed")
