#!/usr/bin/env Rscript

source("oligo_designer.R")

assert_true <- function(value, message) {
  if (!isTRUE(value)) {
    stop(message, call. = FALSE)
  }
}

assert_error <- function(expression, pattern) {
  error <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = identity
  )
  assert_true(!is.null(error), paste("Expected error:", pattern))
  assert_true(
    grepl(pattern, conditionMessage(error)),
    sprintf(
      "Error '%s' does not match '%s'",
      conditionMessage(error),
      pattern
    )
  )
}

base_args <- c(
  "--genome", "genome.fasta",
  "--genome-annotation", "genome.gff",
  "--target-plasmid", "target.fasta",
  "--output-dir", "output",
  "--cds", "geneA"
)

defaults <- parse_designer_args(base_args)
assert_true(
  defaults$annotation_format == "gff",
  "annotation default is not gff"
)
assert_true(
  defaults$chopchop_script == "chopchop/chopchop.py",
  "CHOPCHOP path default is incorrect"
)
assert_true(
  defaults$primer3 == "primer3/src/primer3_core",
  "Primer3 path default is incorrect"
)
assert_true(
  defaults$virtual_pcr_jar == "virtualPCR/dist/virtualPCR.jar",
  "virtualPCR path default is incorrect"
)
assert_true(defaults$n20_mn == 1L, "n20_mn default is not 1")
assert_true(defaults$n20_strands == "random", "strand default is not random")
assert_true(identical(defaults$n20_offtarget, 0L), "MM default is not 0")
assert_true(!defaults$cds_fs && !defaults$ncrna_fs, "FS defaults are not FALSE")
assert_true(
  identical(unname(defaults$left_arm), c(300L, 350L, 400L)),
  "left arm defaults are incorrect"
)
assert_true(
  identical(unname(defaults$right_arm), c(400L, 450L, 500L)),
  "right arm defaults are incorrect"
)
assert_true(
  defaults$n20_arm_min_distance == 40L,
  "N20-to-arm distance default is not 40"
)

configured <- parse_designer_args(c(
  base_args,
  "--n20-mn", "3",
  "--n20-strands", "both",
  "--n20-offtarget", "0,2,4",
  "--cds-fs"
))
assert_true(configured$n20_mn == 3L, "n20_mn was not parsed")
assert_true(configured$n20_strands == "both", "strand mode was not parsed")
assert_true(
  identical(configured$n20_offtarget, c(0L, 2L, 4L)),
  "MM thresholds were not parsed"
)
assert_true(configured$cds_fs, "cds_fs flag was not parsed")

duplicate_targets <- base_args
duplicate_targets[[length(duplicate_targets)]] <- "geneA,geneA"
assert_true(
  identical(parse_designer_args(duplicate_targets)$cds, "geneA"),
  "Duplicate targets were not removed"
)

assert_error(
  parse_designer_args(c(base_args, "--n20-offtarget", "0,-1")),
  "неотрицательные"
)
assert_error(
  parse_designer_args(c(base_args, "--left-arm-min", "401")),
  "min <= opt <= max"
)

feature <- list(start = 100L, end = 299L, length = 200L)
grna_path <- tempfile(fileext = ".tsv")
write_tsv(
  data.frame(
    Rank = 1:4,
    `Target sequence` = rep("AAAAAAAAAAAAAAAAAAAATGG", 4),
    `Genomic location` = paste0("contig:", c(150, 160, 170, 180)),
    Strand = c("+", "-", "+", "-"),
    `Self-complementarity` = c(0, 0, 0, 1),
    MM0 = c(0, 0, 1, 0),
    MM1 = c(1, 3, 0, 0),
    check.names = FALSE
  ),
  grna_path
)
filtered <- filter_grnas(grna_path, feature, "ncrna", c(0L, 2L))
assert_true(nrow(filtered) == 1L, "MM/self-complementarity filtering failed")
assert_true(filtered$strand[[1]] == "+", "Unexpected N20 survived filtering")
assert_error(
  filter_grnas(grna_path, feature, "ncrna", c(0L, 1L, 2L)),
  "только 2 колонок MM"
)

pool <- data.frame(
  rank = 1:4,
  strand = c("+", "-", "-", "+"),
  genomic_location = paste0("contig:", 1:4),
  n20_start = 1:4,
  n20_end = 20:23
)
original_pool <- pool
visited <- 0L
fallback <- visit_grna_sets(pool, 2L, "both", function(candidate) {
  visited <<- visited + 1L
  if (identical(candidate$rank, c(1L, 3L))) candidate else NULL
})
assert_true(visited == 2L, "N20 fallback did not visit the next valid set")
assert_true(
  identical(fallback$rank, c(1L, 3L)),
  "N20 fallback returned the wrong set"
)
assert_true(identical(pool, original_pool), "N20 pool was mutated")

left <- data.frame(
  genome_start = 50L,
  genome_end = 100L,
  arm = "left"
)
right <- data.frame(
  genome_start = 201L,
  genome_end = 250L,
  arm = "right"
)
selected <- list(n20_range = c(141L, 160L))
pair <- choose_pair(
  list(left = left, right = right),
  TRUE,
  list(length = 500L),
  selected,
  minimum_n20_distance = 40L,
  restrict_frame_shift = FALSE
)
assert_true(!is.null(pair), "Exact 40 nt N20 distance was rejected")
pair_with_frame_limit <- choose_pair(
  list(left = left, right = right),
  TRUE,
  list(length = 500L),
  selected,
  minimum_n20_distance = 40L,
  restrict_frame_shift = TRUE
)
assert_true(
  is.null(pair_with_frame_limit),
  "Non-divisible deletion passed frame restriction"
)
minus_pair <- choose_pair(
  list(
    left = data.frame(genome_start = 201L, genome_end = 250L),
    right = data.frame(genome_start = 50L, genome_end = 100L)
  ),
  FALSE,
  list(length = 500L),
  selected,
  minimum_n20_distance = 40L,
  restrict_frame_shift = FALSE
)
assert_true(!is.null(minus_pair), "Minus-strand N20 distance was rejected")

layout <- output_layout("results")
assert_true(
  identical(layout$wet_lab, file.path("results", "WetLab")),
  "WetLab root is incorrect"
)
assert_true(
  identical(layout$tech_report, file.path("results", "TechReport")),
  "TechReport root is incorrect"
)

buffer <- primer3_buffer_parameters()
assert_true(
  identical(
    unname(buffer),
    c(50, 1.5, 0.6, 50)
  ),
  "Primer3 buffer parameters are incorrect"
)
settings_path <- tempfile("primer3-settings-", fileext = ".txt")
write_primer3_settings(settings_path, 400L, 500L, buffer)
settings_text <- readLines(settings_path)
assert_true(
  any(settings_text == "PRIMER_SALT_MONOVALENT=50"),
  "Monovalent salt concentration is absent from Primer3 settings"
)
assert_true(
  any(settings_text == "PRIMER_SALT_DIVALENT=1.5"),
  "Divalent salt concentration is absent from Primer3 settings"
)
assert_true(
  any(settings_text == "PRIMER_DNTP_CONC=0.6"),
  "dNTP concentration is absent from Primer3 settings"
)
assert_true(
  any(settings_text == "PRIMER_DNA_CONC=50"),
  "DNA concentration is absent from Primer3 settings"
)

screening_sizes <- calculate_screening_product_sizes(
  data.frame(PRIMER_PAIR_PRODUCT_SIZE = 500L),
  DNAString(paste(rep("A", 1000L), collapse = "")),
  DNAStringSet(c(
    edited_genome = paste(rep("A", 850L), collapse = "")
  ))
)
assert_true(
  identical(
    screening_sizes,
    c(unsuccessful_insertion_bp = 500L, successful_insertion_bp = 350L)
  ),
  "Screening PCR sizes do not reflect the edited-genome length change"
)

wet_lab_dir <- tempfile("2pac-wet-lab-")
wet_lab_sequences <- DNAStringSet(c(
  geneA_LF = "ACGTACGT",
  geneA_scrF = "TGCATGCA"
))
write_wet_lab_outputs(
  wet_lab_dir,
  list(query_name = "geneA"),
  "cds",
  wet_lab_sequences,
  c("left_arm_forward_primer", "screening_forward_primer"),
  data.frame(
    name = names(wet_lab_sequences),
    purpose = c("left_arm_forward_primer", "screening_forward_primer"),
    annealing_sequence = as.character(wet_lab_sequences),
    tm_c = c(61.2, 63.4),
    stringsAsFactors = FALSE
  ),
  screening_sizes
)
assert_true(
  all(c(
    "final_sequences.fasta",
    "final_sequences.txt",
    "wet_lab_report.txt"
  ) %in% list.files(wet_lab_dir)),
  "WetLab lacks a required file"
)
assert_true(
  !any(c("design.log", "n20_table.tsv", "report.tsv") %in%
    list.files(wet_lab_dir)),
  "WetLab contains technical pipeline output"
)
assert_true(
  length(readDNAStringSet(file.path(wet_lab_dir, "final_sequences.fasta"))) == 2L,
  "WetLab FASTA does not contain the complete final sequence set"
)
wet_lab_table <- read_tsv(
  file.path(wet_lab_dir, "final_sequences.txt"),
  show_col_types = FALSE
)
assert_true(
  identical(wet_lab_table$name, names(wet_lab_sequences)),
  "WetLab TXT does not contain the complete final sequence set"
)
assert_true(
  identical(wet_lab_table$sequence, as.character(wet_lab_sequences)),
  "WetLab FASTA and TXT contain different sequences"
)
wet_lab_report <- readLines(
  file.path(wet_lab_dir, "wet_lab_report.txt"),
  encoding = "UTF-8"
)
assert_true(
  any(grepl("Неуспешная вставка.*500", wet_lab_report)),
  "WetLab report lacks the unsuccessful-insertion PCR size"
)
assert_true(
  any(grepl("Успешная вставка.*350", wet_lab_report)),
  "WetLab report lacks the successful-insertion PCR size"
)
assert_true(
  any(grepl("geneA_LF.*61.2", wet_lab_report)),
  "WetLab report lacks a primer annealing temperature"
)

error_root <- tempfile("2pac-error-")
dir.create(error_root)
stale_wet_lab_dir <- file.path(
  error_root,
  "WetLab",
  "missing_gene_results"
)
dir.create(stale_wet_lab_dir, recursive = TRUE)
stale_wet_lab_files <- file.path(
  stale_wet_lab_dir,
  c("final_sequences.fasta", "final_sequences.txt", "wet_lab_report.txt")
)
file.create(stale_wet_lab_files)
error_input <- list(
  output_dir = error_root,
  annotation = data.frame(
    locus_tag = character(),
    gene = character()
  )
)
target_error <- tryCatch(
  {
    design_target(error_input, "genome", "missing_gene", "cds")
    NULL
  },
  error = identity
)
assert_true(
  inherits(target_error, "target_design_error"),
  "Target error was not classified"
)
target_error_dir <- file.path(
  error_root,
  "TechReport",
  "missing_gene_results"
)
assert_true(
  file.exists(file.path(target_error_dir, "design.log")),
  "design.log was not created"
)
assert_true(
  file.exists(file.path(target_error_dir, "error.txt")),
  "error.txt was not created"
)
error_text <- readLines(file.path(target_error_dir, "error.txt"))
assert_true(any(grepl("stage\\tfeature_lookup", error_text)), "Stage is absent")
assert_true(any(grepl("gene\\tmissing_gene", error_text)), "Gene is absent")
assert_true(
  !any(file.exists(stale_wet_lab_files)),
  "A failed target left stale WetLab output"
)

message("Unit tests passed")
