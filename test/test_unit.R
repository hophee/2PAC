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

error_root <- tempfile("2pac-error-")
dir.create(error_root)
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
target_error_dir <- file.path(error_root, "missing_gene_results")
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

message("Unit tests passed")
