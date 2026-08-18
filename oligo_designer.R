#!/usr/bin/env Rscript
# Batch oligo designer. Functions in this file can also be imported by Shiny.

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(argparser))

source("callPrimer3.R")

extract_gff_attribute <- function(attributes, attribute_name) {
  pattern <- paste0("(?:^|;)\\s*", attribute_name, "=([^;]+)")
  matches <- regexec(pattern, attributes, perl = TRUE)
  values <- regmatches(attributes, matches)
  vapply(
    values,
    function(x) {
      if (length(x) < 2) NA_character_ else utils::URLdecode(trimws(x[[2]]))
    },
    character(1)
  )
}

read_genome_annotation <- function(path, format = "bakta") {
  format <- tolower(format)
  if (format == "bakta") {
    lines <- readLines(path)
    hash_lines <- grep("^#", lines)
    skip_lines <- if (length(hash_lines)) max(hash_lines) - 1 else 0
    annotation <- read_tsv(path, skip = skip_lines, show_col_types = FALSE) |>
      janitor::clean_names()
  } else if (format == "gff") {
    gff <- ape::read.gff(path)
    required <- c("seqid", "type", "start", "end", "strand", "attributes")
    missing <- setdiff(required, names(gff))
    if (length(missing)) {
      stop(
        sprintf(
          "GFF не содержит обязательные колонки: %s",
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    annotation <- data.frame(
      seqid = as.character(gff$seqid),
      type = as.character(gff$type),
      start = as.numeric(gff$start),
      stop = as.numeric(gff$end),
      strand = as.character(gff$strand),
      gene = extract_gff_attribute(gff$attributes, "gene"),
      locus_tag = extract_gff_attribute(gff$attributes, "locus_tag"),
      stringsAsFactors = FALSE
    )
  } else {
    stop(
      "Неподдерживаемый формат аннотации. Допустимы: tsv (bakta), gff",
      call. = FALSE
    )
  }
  required <- c("type", "start", "stop", "strand", "gene", "locus_tag")
  missing <- setdiff(required, names(annotation))
  if (length(missing)) {
    stop(
      sprintf(
        "Аннотация не содержит обязательные колонки: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!"seqid" %in% names(annotation) && "sequence_id" %in% names(annotation)) {
    annotation$seqid <- as.character(annotation$sequence_id)
  }
  if (!"seqid" %in% names(annotation)) {
    stop(
      "Аннотация не содержит колонку seqid или sequence_id",
      call. = FALSE
    )
  }
  annotation
}

find_target_feature <- function(annotation, target_name) {
  by_locus_tag <- which(
    !is.na(annotation$locus_tag) & annotation$locus_tag == target_name
  )
  if (length(by_locus_tag)) {
    return(by_locus_tag[[1]])
  }
  by_gene <- which(!is.na(annotation$gene) & annotation$gene == target_name)
  if (length(by_gene)) {
    return(by_gene[[1]])
  }
  stop(
    sprintf("Не найден ген/feature с locus_tag или gene: %s", target_name),
    call. = FALSE
  )
}

parse_designer_args <- function(args) {
  parser <- arg_parser("Batch oligo designer")
  parser <- add_argument(
    parser,
    "--genome",
    help = "Path to the genome FASTA file"
  )
  parser <- add_argument(
    parser,
    "--genome-annotation",
    help = "Path to the genome annotation file"
  )
  parser <- add_argument(
    parser,
    "--target-plasmid",
    help = "Path to the target plasmid FASTA file"
  )
  parser <- add_argument(
    parser,
    "--output-dir",
    help = "Directory for design results"
  )
  parser <- add_argument(
    parser,
    "--cds",
    help = "CDS targets (comma-separated or space-separated)",
    nargs = Inf
  )
  parser <- add_argument(
    parser,
    "--ncrna",
    help = "ncRNA targets (comma-separated or space-separated)",
    nargs = Inf
  )
  parser <- add_argument(
    parser,
    "--cas-plasmid",
    help = "Path to the Cas plasmid FASTA file"
  )
  parser <- add_argument(
    parser,
    "--annotation-format",
    help = "Annotation format: bakta or gff",
    default = "gff"
  )
  parser <- add_argument(
    parser,
    "--chopchop-script",
    help = "Path to chopchop.py",
    default = "chopchop/chopchop.py"
  )
  parser <- add_argument(
    parser,
    "--primer3",
    help = "Path to primer3_core",
    default = "primer3/src/primer3_core"
  )
  parser <- add_argument(
    parser,
    "--virtual-pcr-jar",
    help = "Path to virtualPCR.jar",
    default = "virtualPCR/dist/virtualPCR.jar"
  )
  parser <- add_argument(
    parser,
    "--n20-mn",
    help = "Minimum number of N20 sequences",
    default = 1L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--n20-strands",
    help = "N20 strand mode: plus, minus, both, or random",
    default = "random"
  )
  parser <- add_argument(
    parser,
    "--n20-offtarget",
    help = "Comma-separated maximum values for MM0, MM1, ...",
    default = "0"
  )
  parser <- add_argument(
    parser,
    "--cds-fs",
    help = "Require the deleted CDS fragment length to be divisible by three",
    flag = TRUE
  )
  parser <- add_argument(
    parser,
    "--ncrna-fs",
    help = "Require the deleted ncRNA fragment length to be divisible by three",
    flag = TRUE
  )
  parser <- add_argument(
    parser,
    "--left-arm-min",
    help = "Minimum left homology arm length",
    default = 300L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--left-arm-opt",
    help = "Preferred initial left homology arm length",
    default = 350L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--left-arm-max",
    help = "Maximum left homology arm length",
    default = 400L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--right-arm-min",
    help = "Minimum right homology arm length",
    default = 400L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--right-arm-opt",
    help = "Preferred initial right homology arm length",
    default = 450L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--right-arm-max",
    help = "Maximum right homology arm length",
    default = 500L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--n20-arm-min-distance",
    help = "Minimum distance from every N20 to both homology arms",
    default = 40L,
    type = "integer"
  )

  normalize_scalar <- function(value) {
    if (is.null(value) || length(value) != 1L || is.na(value[[1]])) {
      return(character())
    }
    value <- trimws(as.character(value[[1]]))
    if (nzchar(value)) value else character()
  }
  normalize_targets <- function(value) {
    if (is.null(value) || !length(value)) {
      return(character())
    }
    value <- as.character(value[!is.na(value)])
    targets <- trimws(unlist(
      strsplit(value, ",", fixed = TRUE),
      use.names = FALSE
    ))
    unique(targets[nzchar(targets)])
  }
  parse_offtarget_thresholds <- function(value) {
    value <- normalize_scalar(value)
    if (!length(value)) {
      stop("--n20-offtarget не может быть пустым", call. = FALSE)
    }
    fields <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
    if (!length(fields) || any(!nzchar(fields))) {
      stop(
        "--n20-offtarget должен быть списком целых чисел через запятую",
        call. = FALSE
      )
    }
    thresholds <- suppressWarnings(as.numeric(fields))
    if (
      anyNA(thresholds) ||
        any(thresholds < 0) ||
        any(thresholds != as.integer(thresholds))
    ) {
      stop(
        "--n20-offtarget допускает только неотрицательные целые числа",
        call. = FALSE
      )
    }
    as.integer(thresholds)
  }
  validate_arm_lengths <- function(side, minimum, preferred, maximum) {
    values <- c(minimum, preferred, maximum)
    if (anyNA(values) || any(values < 1L) || !identical(values, sort(values))) {
      stop(
        sprintf(
          "Для %s плеча требуется 0 < min <= opt <= max",
          side
        ),
        call. = FALSE
      )
    }
    as.integer(values)
  }

  parsed <- parse_args(parser, args)
  annotation_format <- tolower(normalize_scalar(parsed$annotation_format))
  if (!annotation_format %in% c("bakta", "gff")) {
    stop(
      "Некорректный --annotation-format. Допустимы: bakta, gff",
      call. = FALSE
    )
  }
  n20_mn <- as.integer(parsed$n20_mn)
  if (length(n20_mn) != 1L || is.na(n20_mn) || n20_mn < 1L) {
    stop("--n20-mn должен быть положительным целым числом", call. = FALSE)
  }
  n20_strands <- tolower(normalize_scalar(parsed$n20_strands))
  strand_aliases <- c(
    "+" = "plus",
    "plus" = "plus",
    "-" = "minus",
    "minus" = "minus",
    "both" = "both",
    "any" = "random",
    "random" = "random"
  )
  if (
    length(n20_strands) != 1L ||
      !n20_strands %in% names(strand_aliases)
  ) {
    stop(
      "Некорректный --n20-strands. Допустимы: plus, minus, both, random",
      call. = FALSE
    )
  }
  n20_strands <- unname(strand_aliases[[n20_strands]])
  left_arm <- validate_arm_lengths(
    "левого",
    parsed$left_arm_min,
    parsed$left_arm_opt,
    parsed$left_arm_max
  )
  right_arm <- validate_arm_lengths(
    "правого",
    parsed$right_arm_min,
    parsed$right_arm_opt,
    parsed$right_arm_max
  )
  n20_arm_min_distance <- as.integer(parsed$n20_arm_min_distance)
  if (
    length(n20_arm_min_distance) != 1L ||
      is.na(n20_arm_min_distance) ||
      n20_arm_min_distance < 0L
  ) {
    stop(
      "--n20-arm-min-distance должен быть неотрицательным целым числом",
      call. = FALSE
    )
  }

  values <- list(
    genome = normalize_scalar(parsed$genome),
    genome_annotation = normalize_scalar(parsed$genome_annotation),
    target_plasmid = normalize_scalar(parsed$target_plasmid),
    output_dir = normalize_scalar(parsed$output_dir),
    cds = normalize_targets(parsed$cds),
    ncrna = normalize_targets(parsed$ncrna),
    cas_plasmid = normalize_scalar(parsed$cas_plasmid),
    annotation_format = annotation_format,
    chopchop_script = normalize_scalar(parsed$chopchop_script),
    primer3 = normalize_scalar(parsed$primer3),
    virtual_pcr_jar = normalize_scalar(parsed$virtual_pcr_jar),
    n20_mn = n20_mn,
    n20_strands = n20_strands,
    n20_offtarget = parse_offtarget_thresholds(parsed$n20_offtarget),
    cds_fs = isTRUE(parsed$cds_fs),
    ncrna_fs = isTRUE(parsed$ncrna_fs),
    left_arm = setNames(left_arm, c("min", "opt", "max")),
    right_arm = setNames(right_arm, c("min", "opt", "max")),
    n20_arm_min_distance = n20_arm_min_distance
  )

  required <- c("genome", "genome_annotation", "target_plasmid", "output_dir")
  required_options <- c(
    "--genome",
    "--genome-annotation",
    "--target-plasmid",
    "--output-dir"
  )
  missing <- required_options[
    vapply(values[required], length, integer(1)) != 1L
  ]
  if (length(missing)) {
    stop(
      sprintf("Не заданы: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!length(values$cds) && !length(values$ncrna)) {
    stop("Нужен хотя бы один из аргументов --cds или --ncrna", call. = FALSE)
  }
  values
}

run_tool <- function(command, args, stdout = "", stderr = "") {
  status <- system2(command, args = args, stdout = stdout, stderr = stderr)
  if (status != 0) {
    stop(
      sprintf("External command failed (%s): %s", status, command),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

primer3_buffer_parameters <- function() {
  c(
    monovalent_salt_mm = 50,
    divalent_salt_mm = 1.5,
    dntp_mm = 0.6,
    dna_nm = 50
  )
}

output_layout <- function(output_dir) {
  list(
    wet_lab = file.path(output_dir, "WetLab"),
    tech_report = file.path(output_dir, "TechReport")
  )
}

write_run_parameters <- function(input, targets, path) {
  parameters <- c(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    genome_file = input$genome_path,
    genome_annotation_file = input$annotation_path,
    annotation_format = input$annotation_format,
    target_plasmid_file = input$target_plasmid,
    cas_plasmid_file = if (is.null(input$cas_plasmid)) NA_character_ else {
      input$cas_plasmid
    },
    cds_targets = paste(targets$gene[targets$class == "cds"], collapse = ","),
    ncrna_targets = paste(
      targets$gene[targets$class == "ncrna"],
      collapse = ","
    ),
    output_directory = input$output_dir,
    genome_index_directory = file.path(
      output_layout(input$output_dir)$tech_report,
      "genome_indexes"
    ),
    chopchop_script = input$tools$chopchop_script,
    chopchop_config_snapshot = file.path(
      output_layout(input$output_dir)$tech_report,
      "chopchop_config.json"
    ),
    primer3_executable = input$tools$primer3,
    primer3_thermodynamic_parameters = input$tools$primer3_config,
    virtual_pcr_jar = input$tools$virtual_pcr_jar,
    n20_count = input$parameters$n20_mn,
    n20_strands = input$parameters$n20_strands,
    n20_offtarget_thresholds = paste(
      input$parameters$n20_offtarget,
      collapse = ","
    ),
    cds_frame_restriction = input$parameters$cds_fs,
    ncrna_frame_restriction = input$parameters$ncrna_fs,
    setNames(
      input$parameters$left_arm,
      paste0("left_arm_", names(input$parameters$left_arm), "_nt")
    ),
    setNames(
      input$parameters$right_arm,
      paste0("right_arm_", names(input$parameters$right_arm), "_nt")
    ),
    n20_arm_min_distance_nt = input$parameters$n20_arm_min_distance,
    setNames(
      input$parameters$primer3_buffer,
      paste0("primer3_buffer_", names(input$parameters$primer3_buffer))
    )
  )
  write_tsv(
    data.frame(
      parameter = names(parameters),
      value = unname(as.character(parameters)),
      stringsAsFactors = FALSE
    ),
    path
  )
  invisible(path)
}

make_design_input <- function(cli) {
  # TODO: only a complete single-contig bacterial genome is supported for now.
  # Multi-contig FASTA/GFF input needs contig-aware sequence extraction.
  genome_set <- readDNAStringSet(cli$genome[[1]], nrec = 1)
  list(
    genome_path = cli$genome[[1]],
    annotation_path = cli$genome_annotation[[1]],
    annotation_format = cli$annotation_format[[1]],
    genome = genome_set[[1]],
    genome_contig = names(genome_set)[[1]],
    annotation = read_genome_annotation(
      cli$genome_annotation[[1]],
      cli$annotation_format[[1]]
    ),
    target_plasmid = cli$target_plasmid[[1]],
    cas_plasmid = if (length(cli$cas_plasmid)) cli$cas_plasmid[[1]] else NULL,
    output_dir = cli$output_dir[[1]],
    tools = list(
      chopchop_script = cli$chopchop_script[[1]],
      primer3 = cli$primer3[[1]],
      primer3_config = "primer3/src/primer3_config",
      virtual_pcr_jar = cli$virtual_pcr_jar[[1]]
    ),
    parameters = list(
      n20_mn = cli$n20_mn,
      n20_strands = cli$n20_strands,
      n20_offtarget = cli$n20_offtarget,
      cds_fs = cli$cds_fs,
      ncrna_fs = cli$ncrna_fs,
      left_arm = cli$left_arm,
      right_arm = cli$right_arm,
      n20_arm_min_distance = cli$n20_arm_min_distance,
      primer3_buffer = primer3_buffer_parameters()
    )
  )
}

feature_record <- function(input, gene_name) {
  idx <- find_target_feature(input$annotation, gene_name)
  row <- input$annotation[idx, , drop = FALSE]
  contig <- as.character(row$seqid[[1]])
  candidates <- input$annotation
  if (!is.na(contig)) {
    candidates <- filter(candidates, seqid == contig)
  }
  previous_stop <- suppressWarnings(max(
    candidates$stop[candidates$stop < row$start[[1]]],
    na.rm = TRUE
  ))
  next_start <- suppressWarnings(min(
    candidates$start[candidates$start > row$stop[[1]]],
    na.rm = TRUE
  ))
  if (!is.finite(previous_stop)) {
    previous_stop <- 0L
  }
  if (!is.finite(next_start)) {
    next_start <- length(input$genome) + 1L
  }
  list(
    index = idx,
    query_name = gene_name,
    display_name = if (!is.na(row$gene[[1]]) && nzchar(row$gene[[1]])) {
      row$gene[[1]]
    } else {
      gene_name
    },
    contig = if (!is.na(contig)) contig else input$genome_contig,
    start = as.integer(row$start[[1]]),
    end = as.integer(row$stop[[1]]),
    strand = as.character(row$strand[[1]]),
    length = as.integer(row$stop[[1]] - row$start[[1]] + 1),
    left_bound = as.integer(previous_stop + 1L),
    right_bound = as.integer(next_start - 1L)
  )
}

cut_interval <- function(feature, design_class, genome_length) {
  if (feature$length <= 500) {
    if (design_class == "ncrna") {
      extension <- round(feature$length / 10, -1) * 2
      interval <- c(feature$start - extension, feature$end + extension)
    } else {
      interval <- c(feature$start, feature$end)
    }
  } else if (feature$length < 1500) {
    midpoint <- feature$start + feature$length %/% 2
    interval <- c(midpoint - 250, midpoint + 250)
  } else {
    flank <- floor((feature$end - feature$start) / 3)
    interval <- c(feature$start + flank, feature$end - flank)
  }
  interval <- pmax(1L, pmin(as.integer(interval), genome_length))
  if (design_class == "ncrna") {
    interval[[1]] <- max(interval[[1]], feature$left_bound)
    interval[[2]] <- min(interval[[2]], feature$right_bound)
  }
  interval
}

prepare_chopchop_assets <- function(input) {
  genome_name <- tools::file_path_sans_ext(basename(input$genome_path))
  index_dir <- file.path(
    output_layout(input$output_dir)$tech_report,
    "genome_indexes"
  )
  dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
  index_dir <- normalizePath(index_dir)
  two_bit <- file.path(index_dir, paste0(genome_name, ".2bit"))
  bowtie_prefix <- file.path(index_dir, genome_name)
  if (!file.exists(two_bit)) {
    run_tool("faToTwoBit", c(input$genome_path, two_bit))
  }
  if (!file.exists(paste0(bowtie_prefix, ".1.ebwt"))) {
    run_tool("bowtie-build", c(input$genome_path, bowtie_prefix))
  }
  list(name = genome_name, directory = index_dir)
}

configure_chopchop <- function(input, genome_assets) {
  script_path <- normalizePath(input$tools$chopchop_script)
  chopchop_dir <- dirname(script_path)
  two_bit_to_fa <- Sys.which("twoBitToFa")
  bowtie <- Sys.which("bowtie")
  if (!nzchar(two_bit_to_fa) || !nzchar(bowtie)) {
    stop("twoBitToFa или bowtie не найдены в PATH", call. = FALSE)
  }
  config <- file.path(chopchop_dir, "config_local.json")
  writeLines(
    c(
      "{",
      "  \"PATH\": {",
      paste0("    \"PRIMER3\": \"", normalizePath(input$tools$primer3), "\","),
      paste0("    \"BOWTIE\": \"", bowtie, "\","),
      paste0("    \"TWOBITTOFA\": \"", two_bit_to_fa, "\","),
      paste0("    \"TWOBIT_INDEX_DIR\": \"", genome_assets$directory, "\","),
      paste0("    \"BOWTIE_INDEX_DIR\": \"", genome_assets$directory, "\","),
      paste0("    \"ISOFORMS_INDEX_DIR\": \"", genome_assets$directory, "\","),
      paste0("    \"ISOFORMS_MT_DIR\": \"", genome_assets$directory, "\","),
      paste0("    \"GENE_TABLE_INDEX_DIR\": \"", genome_assets$directory, "\""),
      "  },",
      "  \"THREADS\": 1",
      "}"
    ),
    config
  )
  invisible(config)
}

run_chopchop <- function(
  input,
  genome_name,
  feature,
  design_class,
  target_dir
) {
  interval <- cut_interval(feature, design_class, length(input$genome))
  target <- paste0(feature$contig, ":", interval[[1]], "-", interval[[2]])
  table_path <- file.path(target_dir, "n20_table.tsv")
  run_tool(
    "python",
    c(
      input$tools$chopchop_script,
      "-Target",
      target,
      "-G",
      genome_name,
      "-M",
      "NGG",
      "-T",
      "1",
      "-g",
      "20",
      "-m",
      "3",
      "--padSize",
      "0",
      "-O",
      "50",
      "--scoringMethod",
      "DOENCH_2016",
      "-o",
      target_dir
    ),
    stdout = table_path,
    stderr = file.path(target_dir, "chopchop.stderr.log")
  )
  if (!file.exists(table_path) || file.size(table_path) == 0) {
    stop("ChopChop did not produce n20_table.tsv", call. = FALSE)
  }
  offtarget_files <- list.files(
    target_dir,
    pattern = "\\.offtargets$",
    full.names = TRUE
  )
  if (length(offtarget_files)) {
    offtarget_dir <- file.path(target_dir, "chopchop_offtargets")
    dir.create(offtarget_dir, showWarnings = FALSE)
    file.rename(
      offtarget_files,
      file.path(offtarget_dir, basename(offtarget_files))
    )
  }
  unlink(file.path(
    target_dir,
    c(
      "sequence.fa",
      "gene_file.fa",
      "output.sam",
      "bowtie.err",
      "twoBitToFa.err"
    )
  ))
  interval
}

filter_grnas <- function(
  table_path,
  feature,
  design_class,
  offtarget_thresholds
) {
  grnas <- read_tsv(table_path, show_col_types = FALSE) |>
    janitor::clean_names()
  required <- c(
    "genomic_location",
    "target_sequence",
    "strand",
    "self_complementarity"
  )
  missing <- setdiff(required, names(grnas))
  if (length(missing)) {
    stop(
      sprintf(
        "Таблица CHOPCHOP не содержит обязательные колонки: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  mm_columns <- grep("^mm[0-9]+$", names(grnas), value = TRUE)
  mm_columns <- mm_columns[order(as.integer(sub("^mm", "", mm_columns)))]
  if (length(offtarget_thresholds) > length(mm_columns)) {
    stop(
      sprintf(
        paste(
          "--n20-offtarget содержит %d значений,",
          "но CHOPCHOP вывел только %d колонок MM"
        ),
        length(offtarget_thresholds),
        length(mm_columns)
      ),
      call. = FALSE
    )
  }
  applied_mm_columns <- head(mm_columns, length(offtarget_thresholds))
  keep_offtarget <- rep(TRUE, nrow(grnas))
  for (i in seq_along(applied_mm_columns)) {
    values <- suppressWarnings(as.numeric(grnas[[applied_mm_columns[[i]]]]))
    keep_offtarget <- keep_offtarget &
      !is.na(values) &
      values <= offtarget_thresholds[[i]]
  }
  grnas <- grnas[keep_offtarget, , drop = FALSE] |>
    mutate(
      genomic_start = as.numeric(sub("^.*:", "", genomic_location)),
      genomic_end = genomic_start + 22L,
      n20_start = genomic_start,
      n20_end = genomic_start + 19L,
      mid_closeness = abs(
        (feature$start + feature$end) %/%
          2 -
          (n20_start + n20_end) %/% 2
      ) /
        feature$length
    ) |>
    filter(self_complementarity == 0)
  if (anyNA(grnas$genomic_start)) {
    stop(
      "Не удалось разобрать genomic_location в таблице CHOPCHOP",
      call. = FALSE
    )
  }
  if (design_class == "cds") {
    grnas <- filter(grnas, mid_closeness <= 0.18)
  }
  arrange(grnas, mid_closeness)
}

prepare_grna_pool <- function(grnas, n20_mn, strand_mode) {
  if (n20_mn > 1L && strand_mode == "plus") {
    grnas <- filter(grnas, strand == "+")
  } else if (n20_mn > 1L && strand_mode == "minus") {
    grnas <- filter(grnas, strand == "-")
  }
  if (nrow(grnas) < n20_mn) {
    stop(
      sprintf(
        "Недостаточно подходящих N20: найдено %d, требуется не менее %d",
        nrow(grnas),
        n20_mn
      ),
      call. = FALSE
    )
  }
  if (
    n20_mn > 1L &&
      strand_mode == "both" &&
      (!any(grnas$strand == "+") || !any(grnas$strand == "-"))
  ) {
    stop(
      "Для режима both нужны подходящие N20 на плюс- и минус-цепях",
      call. = FALSE
    )
  }
  grnas
}

visit_grna_sets <- function(grnas, n20_mn, strand_mode, visitor) {
  selected_indices <- integer(n20_mn)
  visit <- function(depth, first_index) {
    if (depth > n20_mn) {
      candidates <- grnas[selected_indices, , drop = FALSE]
      if (
        n20_mn > 1L &&
          strand_mode == "both" &&
          !all(c("+", "-") %in% candidates$strand)
      ) {
        return(NULL)
      }
      return(visitor(candidates))
    }
    remaining <- n20_mn - depth
    last_index <- nrow(grnas) - remaining
    for (index in seq.int(first_index, last_index)) {
      selected_indices[[depth]] <<- index
      result <- visit(depth + 1L, index + 1L)
      if (!is.null(result)) {
        return(result)
      }
    }
    NULL
  }
  visit(1L, 1L)
}

write_selected_grnas <- function(selected, feature, target_dir) {
  write_tsv(
    selected$table,
    file.path(target_dir, "selected_n20_table.tsv")
  )
  n20 <- DNAStringSet(substr(selected$table$target_sequence, 1, 20))
  names(n20) <- paste0(feature$display_name, "_n20_", seq_along(n20))
  writeXStringSet(n20, file.path(target_dir, "selected_n20.fasta"))
  invisible(selected)
}

write_primer3_settings <- function(
  path,
  left_length,
  right_length,
  buffer = primer3_buffer_parameters()
) {
  product <- round(c(
    min(left_length, right_length),
    max(left_length, right_length) * 1.5
  ))
  writeLines(
    c(
      "Primer3 File - http://primer3.org",
      "P3_FILE_TYPE=settings",
      "",
      "PRIMER_TASK=generic",
      "PRIMER_PICK_LEFT_PRIMER=1",
      "PRIMER_PICK_RIGHT_PRIMER=1",
      "PRIMER_NUM_RETURN=10",
      "PRIMER_MIN_SIZE=18",
      "PRIMER_OPT_SIZE=21",
      "PRIMER_MAX_SIZE=27",
      "PRIMER_MIN_TM=59.0",
      "PRIMER_OPT_TM=60.0",
      "PRIMER_MAX_TM=61.0",
      paste0("PRIMER_SALT_MONOVALENT=", buffer[["monovalent_salt_mm"]]),
      paste0("PRIMER_SALT_DIVALENT=", buffer[["divalent_salt_mm"]]),
      paste0("PRIMER_DNTP_CONC=", buffer[["dntp_mm"]]),
      paste0("PRIMER_DNA_CONC=", buffer[["dna_nm"]]),
      "PRIMER_PAIR_MAX_DIFF_TM=8.0",
      "PRIMER_MIN_GC=40.0",
      "PRIMER_MAX_GC=60.0",
      "PRIMER_MAX_SELF_ANY=12.0",
      "PRIMER_MAX_SELF_END=8.0",
      "PRIMER_PAIR_MAX_COMPL_ANY=12.0",
      "PRIMER_PAIR_MAX_COMPL_END=8.0",
      "PRIMER_MAX_POLY_X=5",
      paste0("PRIMER_PRODUCT_SIZE_RANGE=", product[[1]], "-", product[[2]]),
      "PRIMER_EXPLAIN_FLAG=1",
      "PRIMER_FIRST_BASE_INDEX=1",
      "="
    ),
    path
  )
}

add_genome_positions <- function(left, right, arm, plus_strand) {
  if (plus_strand) {
    left <- mutate(
      left,
      genome_start = PRIMER_LEFT_pos + arm$left_start - 1,
      genome_end = PRIMER_RIGHT_pos + arm$left_start - 1
    )
    right <- mutate(
      right,
      genome_start = PRIMER_LEFT_pos + arm$right_start - 1,
      genome_end = PRIMER_RIGHT_pos + arm$right_start - 1
    )
  } else {
    left <- mutate(
      left,
      genome_start = arm$left_end - PRIMER_RIGHT_pos - 1,
      genome_end = arm$left_end - PRIMER_LEFT_pos - 1
    )
    right <- mutate(
      right,
      genome_start = arm$right_end - PRIMER_RIGHT_pos + 1,
      genome_end = arm$right_end - PRIMER_LEFT_pos + 1
    )
  }
  list(left = left, right = right)
}

choose_pair <- function(
  primers,
  plus_strand,
  feature,
  selected,
  minimum_n20_distance,
  restrict_frame_shift
) {
  left <- primers$left
  right <- primers$right
  if (!nrow(left) || !nrow(right)) {
    return(NULL)
  }
  if (plus_strand) {
    deleted <- outer(right$genome_start, left$genome_end, "-") - 1
    lower_boundary <- matrix(
      rep(left$genome_end, each = nrow(right)),
      nrow = nrow(right)
    )
    upper_boundary <- matrix(
      rep(right$genome_start, nrow(left)),
      nrow = nrow(right)
    )
  } else {
    deleted <- outer(left$genome_start, right$genome_end, "-") - 1
    lower_boundary <- matrix(
      rep(right$genome_end, each = nrow(left)),
      nrow = nrow(left)
    )
    upper_boundary <- matrix(
      rep(left$genome_start, nrow(right)),
      nrow = nrow(left)
    )
  }
  left_distance <- selected$n20_range[[1]] - lower_boundary - 1L
  right_distance <- upper_boundary - selected$n20_range[[2]] - 1L
  eligible <- deleted >= 0L &
    deleted < feature$length &
    left_distance >= minimum_n20_distance &
    right_distance >= minimum_n20_distance
  if (restrict_frame_shift) {
    eligible <- eligible & deleted %% 3 == 0
  }
  indices <- which(eligible)
  if (!length(indices)) {
    return(NULL)
  }
  index <- indices[[which.min(deleted[indices])]]
  rc <- arrayInd(index, dim(deleted))
  if (plus_strand) {
    bind_rows(left[rc[[2]], ], right[rc[[1]], ])
  } else {
    bind_rows(left[rc[[1]], ], right[rc[[2]], ])
  }
}

design_homology_arms <- function(
  input,
  feature,
  selected,
  design_class,
  target_dir
) {
  interval <- cut_interval(feature, design_class, length(input$genome))
  parameters <- input$parameters
  left_limits <- parameters$left_arm
  right_limits <- parameters$right_arm
  minimum_n20_distance <- parameters$n20_arm_min_distance
  restrict_frame_shift <- if (design_class == "cds") {
    parameters$cds_fs
  } else {
    parameters$ncrna_fs
  }
  settings <- file.path(target_dir, "primer3_settings.txt")
  write_primer3_settings(
    settings,
    left_limits[["max"]],
    right_limits[["max"]],
    input$parameters$primer3_buffer
  )
  if (!file.exists(input$tools$primer3)) {
    stop("primer3_core не найден", call. = FALSE)
  }

  for (attempt in seq_len(1000)) {
    left_length <- min(
      left_limits[["opt"]] + attempt - 1L,
      left_limits[["max"]]
    )
    right_length <- min(
      right_limits[["opt"]] + attempt - 1L,
      right_limits[["max"]]
    )
    if (feature$strand == "+") {
      left_end <- interval[[1]] - 1L
      left_start <- max(1L, left_end - left_length + 1L)
      right_start <- interval[[2]] + 1L
      right_end <- min(
        length(input$genome),
        right_start + right_length - 1L
      )
      left_seq <- input$genome[left_start:left_end]
      right_seq <- input$genome[right_start:right_end]
    } else {
      right_end <- interval[[1]] - 1L
      left_start <- interval[[2]] + 1L
      left_end <- min(
        length(input$genome),
        left_start + left_length - 1L
      )
      right_start <- max(1L, right_end - right_length + 1L)
      left_seq <- complement(reverse(input$genome[left_start:left_end]))
      right_seq <- complement(reverse(input$genome[right_start:right_end]))
    }
    if (
      length(left_seq) < left_limits[["min"]] ||
        length(right_seq) < right_limits[["min"]]
    ) {
      stop(
        "Граница генома не позволяет выдержать минимальную длину плеч",
        call. = FALSE
      )
    }
    writeXStringSet(
      DNAStringSet(c(left_arm = left_seq, right_arm = right_seq)),
      file.path(target_dir, "homology_arms_before_primer_search.fasta")
    )
    tm <- if (design_class == "cds") c(62.5, 63, 63.5) else c(60.5, 61, 62.5)
    left <- tryCatch(
      callPrimer3(
        as.character(left_seq),
        paste0(left_limits[["min"]], "-", length(left_seq)),
        tm,
        2,
        "left_arm",
        primer_num = 10,
        primer3 = input$tools$primer3,
        thermo.param = input$tools$primer3_config,
        settings = settings,
        report = file.path(target_dir, "left_arm_report.txt")
      ),
      error = function(e) NULL
    )
    right <- tryCatch(
      callPrimer3(
        as.character(right_seq),
        paste0(right_limits[["min"]], "-", length(right_seq)),
        tm,
        2,
        "right_arm",
        primer_num = 10,
        primer3 = input$tools$primer3,
        thermo.param = input$tools$primer3_config,
        settings = settings,
        report = file.path(target_dir, "right_arm_report.txt")
      ),
      error = function(e) NULL
    )
    if (is.data.frame(left) && is.data.frame(right)) {
      arm <- list(
        left_start = left_start,
        left_end = left_end,
        right_start = right_start,
        right_end = right_end
      )
      positions <- add_genome_positions(left, right, arm, feature$strand == "+")
      pair <- choose_pair(
        positions,
        feature$strand == "+",
        feature,
        selected,
        minimum_n20_distance,
        restrict_frame_shift
      )
      if (!is.null(pair)) {
        pair_arm_lengths <- pair$PRIMER_RIGHT_pos -
          pair$PRIMER_LEFT_pos +
          1L
        valid_arm_lengths <- pair_arm_lengths[[1]] >= left_limits[["min"]] &&
          pair_arm_lengths[[1]] <= left_limits[["max"]] &&
          pair_arm_lengths[[2]] >= right_limits[["min"]] &&
          pair_arm_lengths[[2]] <= right_limits[["max"]]
        ticks <- sort(unlist(pair[, c("genome_start", "genome_end")]))
        if (
          valid_arm_lengths &&
            ticks[[2]] >= feature$start &&
            ticks[[3]] <= feature$end
        ) {
          write_tsv(
            bind_rows(positions$left, positions$right),
            file.path(target_dir, "primer3_table.tsv")
          )
          return(list(
            pair = pair,
            left = left_seq,
            right = right_seq,
            arm = arm,
            ticks = ticks
          ))
        }
      }
    }
    interval <- interval + c(1L, -1L)
    required_deletion <- selected$n20_range +
      c(-minimum_n20_distance - 1L, minimum_n20_distance + 1L)
    if (
      interval[[1]] > required_deletion[[1]] ||
        interval[[2]] < required_deletion[[2]] ||
        diff(interval) + 1 < feature$length * .3
    ) {
      break
    }
  }
  NULL
}

calculate_screening_product_sizes <- function(
  screening,
  original_genome,
  edited_genome
) {
  reference_size <- if (
    "PRIMER_PAIR_PRODUCT_SIZE" %in% names(screening)
  ) {
    suppressWarnings(as.integer(screening$PRIMER_PAIR_PRODUCT_SIZE[[1]]))
  } else {
    as.integer(
      screening$genome_end[[1]] - screening$genome_start[[1]] + 1L
    )
  }
  edited_genome_size <- if (inherits(edited_genome, "XStringSet")) {
    width(edited_genome)[[1]]
  } else {
    length(edited_genome)
  }
  edited_size <- reference_size + edited_genome_size - length(original_genome)
  if (
    length(reference_size) != 1L ||
      is.na(reference_size) ||
      reference_size < 1L ||
      length(edited_size) != 1L ||
      is.na(edited_size) ||
      edited_size < 1L
  ) {
    stop("Не удалось рассчитать размеры скрининговых ПЦР-продуктов", call. = FALSE)
  }
  c(
    unsuccessful_insertion_bp = reference_size,
    successful_insertion_bp = as.integer(edited_size)
  )
}

write_wet_lab_outputs <- function(
  wet_lab_dir,
  feature,
  design_class,
  sequences,
  sequence_purposes,
  primer_metrics,
  screening_product_sizes
) {
  required_metrics <- c("name", "purpose", "annealing_sequence", "tm_c")
  if (
    !length(sequences) ||
      length(sequence_purposes) != length(sequences) ||
      is.null(names(sequences)) ||
      anyNA(names(sequences)) ||
      any(!nzchar(names(sequences))) ||
      !all(required_metrics %in% names(primer_metrics)) ||
      !all(
        c("unsuccessful_insertion_bp", "successful_insertion_bp") %in%
          names(screening_product_sizes)
      )
  ) {
    stop("Неполный набор результатов для WetLab", call. = FALSE)
  }
  dir.create(wet_lab_dir, recursive = TRUE, showWarnings = FALSE)
  writeXStringSet(
    sequences,
    file.path(wet_lab_dir, "final_sequences.fasta")
  )
  sequence_table <- data.frame(
    name = names(sequences),
    sequence = as.character(sequences),
    purpose = unname(sequence_purposes),
    stringsAsFactors = FALSE
  )
  write_tsv(
    sequence_table,
    file.path(wet_lab_dir, "final_sequences.txt")
  )

  sequence_lines <- apply(
    sequence_table,
    1L,
    function(row) paste(row, collapse = "\t")
  )
  tm_table <- primer_metrics
  tm_table$tm_c <- format(
    round(tm_table$tm_c, 1),
    nsmall = 1,
    trim = TRUE
  )
  tm_lines <- apply(
    tm_table,
    1L,
    function(row) paste(row, collapse = "\t")
  )
  report <- c(
    "2PAC: отчёт для мокрой лаборатории",
    paste("Цель", feature$query_name, sep = "\t"),
    paste("Класс", design_class, sep = "\t"),
    "",
    "Итоговый набор последовательностей",
    paste(names(sequence_table), collapse = "\t"),
    sequence_lines,
    "",
    paste(
      "Температуры отжига праймеров",
      "Tm рассчитана для участка отжига без сервисных последовательностей",
      sep = "\n"
    ),
    paste(names(tm_table), collapse = "\t"),
    tm_lines,
    "",
    "Размеры скрининговых ПЦР-продуктов",
    paste(
      "Неуспешная вставка (исходный аллель), п.н.",
      screening_product_sizes[["unsuccessful_insertion_bp"]],
      sep = "\t"
    ),
    paste(
      "Успешная вставка (редактированный аллель), п.н.",
      screening_product_sizes[["successful_insertion_bp"]],
      sep = "\t"
    )
  )
  writeLines(
    enc2utf8(report),
    file.path(wet_lab_dir, "wet_lab_report.txt"),
    useBytes = TRUE
  )
  invisible(wet_lab_dir)
}

write_design_outputs <- function(
  input,
  feature,
  selected,
  arms,
  design_class,
  target_dir
) {
  pair <- arms$pair
  gap <- arms$ticks[[3]] - arms$ticks[[2]] - 1L
  bridge <- "ATGACTGCCCGCAAG"
  restrict_frame_shift <- if (design_class == "cds") {
    input$parameters$cds_fs
  } else {
    input$parameters$ncrna_fs
  }
  frame_status <- if (restrict_frame_shift) {
    "divisible_by_three"
  } else {
    "not_restricted"
  }
  if (design_class == "cds") {
    bridge <- substr(bridge, 1, 15 - (gap %% 3))
    if (!restrict_frame_shift && gap %% 3 != 0) {
      frame_status <- sprintf("bridge shortened to %d nt", nchar(bridge))
    }
  }
  bridge_rc <- as.character(complement(reverse(DNAString(bridge))))
  final_arms <- DNAStringSet(c(
    left_homology_arm = arms$left[
      pair$PRIMER_LEFT_pos[[1]]:pair$PRIMER_RIGHT_pos[[1]]
    ],
    right_homology_arm = arms$right[
      pair$PRIMER_LEFT_pos[[2]]:pair$PRIMER_RIGHT_pos[[2]]
    ]
  ))
  writeXStringSet(final_arms, file.path(target_dir, "homology_arms.fasta"))
  sgrnas <- DNAStringSet(paste0(
    "ACGACTAGT",
    substr(selected$table$target_sequence, 1, 20),
    "GTTTTAGAGCTAGAAATAGCAAGTTaaaataaggct"
  ))
  names(sgrnas) <- paste0(feature$display_name, "_sgF", seq_along(sgrnas))
  arm_primers <- DNAStringSet(c(
    paste0("AGCGTCAACT", pair$PRIMER_LEFT_SEQUENCE[[1]]),
    paste0(bridge_rc, pair$PRIMER_RIGHT_SEQUENCE[[1]]),
    paste0(bridge, pair$PRIMER_LEFT_SEQUENCE[[2]]),
    paste0("ACGCTGCAG", pair$PRIMER_RIGHT_SEQUENCE[[2]])
  ))
  names(arm_primers) <- paste0(
    feature$display_name,
    c("_LF", "_LR", "_RF", "_RR")
  )
  sgrna_reverse <- DNAStringSet("AGTTGACGCTAAAAAAAGCACCGACTCGGTGCC")
  names(sgrna_reverse) <- paste0(feature$display_name, "_sgR")
  all_primers <- c(sgrnas, sgrna_reverse, arm_primers)
  writeXStringSet(all_primers, file.path(target_dir, "all_primers.fasta"))
  plain <- DNAStringSet(c(
    substr(selected$table$target_sequence, 1, 20),
    pair$PRIMER_LEFT_SEQUENCE[[1]],
    pair$PRIMER_RIGHT_SEQUENCE[[1]],
    pair$PRIMER_LEFT_SEQUENCE[[2]],
    pair$PRIMER_RIGHT_SEQUENCE[[2]]
  ))
  names(plain) <- c(names(sgrnas), names(arm_primers))
  plain_path <- file.path(target_dir, "primers_without_service_sequences.fasta")
  writeXStringSet(plain, plain_path)

  offtarget_range <- range(unlist(pair[, c("genome_start", "genome_end")]))
  screening_range <- pmax(
    1L,
    pmin(as.integer(offtarget_range + c(-200L, 200L)), length(input$genome))
  )
  screening_seq <- input$genome[screening_range[[1]]:screening_range[[2]]]
  screening <- tryCatch(
    callPrimer3(
      as.character(screening_seq),
      paste0(length(screening_seq) - 100L, "-", length(screening_seq)),
      c(62.5, 63, 63.5),
      2,
      "genome_screening",
      primer_num = 5,
      primer3 = input$tools$primer3,
      thermo.param = input$tools$primer3_config,
      settings = file.path(target_dir, "primer3_settings.txt"),
      report = file.path(target_dir, "genome_screening_report.txt")
    ),
    error = function(e) NULL
  )
  if (!is.data.frame(screening) || !nrow(screening)) {
    stop("Не удалось подобрать скрининговые праймеры", call. = FALSE)
  }
  screening <- mutate(
    screening,
    genome_start = PRIMER_LEFT_pos + screening_range[[1]] - 1L,
    genome_end = PRIMER_RIGHT_pos + screening_range[[1]] - 1L
  )
  write_tsv(screening, file.path(target_dir, "screening_primer3_table.tsv"))
  screening_primers <- DNAStringSet(c(
    screening$PRIMER_LEFT_SEQUENCE[[1]],
    screening$PRIMER_RIGHT_SEQUENCE[[1]]
  ))
  names(screening_primers) <- paste0(feature$display_name, "_scr", c("F", "R"))
  screening_path <- file.path(target_dir, "screening_primers.fasta")
  writeXStringSet(screening_primers, screening_path)

  starts <- pair$PRIMER_LEFT_pos + pair$PRIMER_LEFT_len
  ends <- pair$PRIMER_RIGHT_pos - pair$PRIMER_RIGHT_len
  prefix <- if (offtarget_range[[1]] > 1L) {
    input$genome[1:offtarget_range[[1]]]
  } else {
    DNAString("")
  }
  suffix <- input$genome[offtarget_range[[2]]:length(input$genome)]
  edited_genome <- DNAStringSet(c(
    edited_genome = c(
      prefix,
      arms$left[starts[[1]]:ends[[1]]],
      arm_primers[[3]],
      arms$right[starts[[2]]:ends[[2]]],
      suffix
    )
  ))
  writeXStringSet(edited_genome, file.path(target_dir, "edited_genome.fasta"))
  screening_product_sizes <- calculate_screening_product_sizes(
    screening,
    input$genome,
    edited_genome
  )
  final_sequences <- c(all_primers, screening_primers)
  primer_purposes <- c(
    "left_arm_forward_primer",
    "left_arm_reverse_primer",
    "right_arm_forward_primer",
    "right_arm_reverse_primer",
    "screening_forward_primer",
    "screening_reverse_primer"
  )
  sequence_purposes <- c(
    rep("sgRNA_forward_oligo", length(sgrnas)),
    "sgRNA_reverse_oligo",
    primer_purposes
  )
  primer_metrics <- data.frame(
    name = c(names(arm_primers), names(screening_primers)),
    purpose = primer_purposes,
    annealing_sequence = c(
      pair$PRIMER_LEFT_SEQUENCE[[1]],
      pair$PRIMER_RIGHT_SEQUENCE[[1]],
      pair$PRIMER_LEFT_SEQUENCE[[2]],
      pair$PRIMER_RIGHT_SEQUENCE[[2]],
      screening$PRIMER_LEFT_SEQUENCE[[1]],
      screening$PRIMER_RIGHT_SEQUENCE[[1]]
    ),
    tm_c = c(
      pair$PRIMER_LEFT_TM[[1]],
      pair$PRIMER_RIGHT_TM[[1]],
      pair$PRIMER_LEFT_TM[[2]],
      pair$PRIMER_RIGHT_TM[[2]],
      screening$PRIMER_LEFT_TM[[1]],
      screening$PRIMER_RIGHT_TM[[1]]
    ),
    stringsAsFactors = FALSE
  )

  writeLines(
    c(
      paste("target", feature$query_name, sep = "\t"),
      paste("class", design_class, sep = "\t"),
      paste("n20_count", nrow(selected$table), sep = "\t"),
      paste(
        "n20_strands",
        paste(sort(unique(selected$table$strand)), collapse = ","),
        sep = "\t"
      ),
      paste(
        "n20_offtarget_thresholds",
        paste(input$parameters$n20_offtarget, collapse = ","),
        sep = "\t"
      ),
      paste(
        "n20_arm_min_distance",
        input$parameters$n20_arm_min_distance,
        sep = "\t"
      ),
      paste("deleted_nt", gap, sep = "\t"),
      paste("frame_status", frame_status, sep = "\t"),
      paste("left_arm_nt", length(final_arms[[1]]), sep = "\t"),
      paste("right_arm_nt", length(final_arms[[2]]), sep = "\t"),
      paste(
        "left_forward_primer_tm",
        round(pair$PRIMER_LEFT_TM[[1]], 1),
        sep = "\t"
      ),
      paste(
        "left_reverse_primer_tm",
        round(pair$PRIMER_RIGHT_TM[[1]], 1),
        sep = "\t"
      ),
      paste(
        "right_forward_primer_tm",
        round(pair$PRIMER_LEFT_TM[[2]], 1),
        sep = "\t"
      ),
      paste(
        "right_reverse_primer_tm",
        round(pair$PRIMER_RIGHT_TM[[2]], 1),
        sep = "\t"
      ),
      paste(
        "screening_forward_tm",
        round(screening$PRIMER_LEFT_TM[[1]], 1),
        sep = "\t"
      ),
      paste(
        "screening_reverse_tm",
        round(screening$PRIMER_RIGHT_TM[[1]], 1),
        sep = "\t"
      ),
      paste(
        "screening_unsuccessful_insertion_bp",
        screening_product_sizes[["unsuccessful_insertion_bp"]],
        sep = "\t"
      ),
      paste(
        "screening_successful_insertion_bp",
        screening_product_sizes[["successful_insertion_bp"]],
        sep = "\t"
      )
    ),
    file.path(target_dir, "report.tsv")
  )
  list(
    all_primers_path = file.path(target_dir, "all_primers.fasta"),
    plain_path = plain_path,
    screening_path = screening_path,
    wet_lab = list(
      sequences = final_sequences,
      sequence_purposes = sequence_purposes,
      primer_metrics = primer_metrics,
      screening_product_sizes = screening_product_sizes
    )
  )
}

run_virtual_pcr <- function(input, target_dir, primer_paths) {
  if (!file.exists(input$tools$virtual_pcr_jar)) {
    stop("virtualPCR JAR не найден", call. = FALSE)
  }
  check <- function(label, targets_path, primers_path, output_name) {
    config <- file.path(target_dir, paste0("virtual_pcr_", label, ".conf"))
    writeLines(
      c(
        paste0("targets_path=", targets_path),
        paste0("output_path=", file.path(target_dir, output_name)),
        paste0("primers_path=", primers_path),
        "type=primer",
        "ShowPCRProducts=true",
        "ShowPrimerAlignment=true",
        "ShowPrimerAlignmentPCRproduct=false",
        "primerstatistic=true"
      ),
      config
    )
    run_tool(
      "java",
      c("-jar", input$tools$virtual_pcr_jar, config),
      stdout = file.path(
        target_dir,
        paste0("virtual_pcr_", label, ".stdout.log")
      ),
      stderr = file.path(
        target_dir,
        paste0("virtual_pcr_", label, ".stderr.log")
      )
    )
    unlink(config)
  }
  check(
    "genome",
    input$genome_path,
    primer_paths$all_primers_path,
    "genome_offtarget_check.txt"
  )
  check(
    "genome_non_service",
    input$genome_path,
    primer_paths$plain_path,
    "genome_offtarget_non_service_seq.txt"
  )
  check(
    "target_plasmid_screening",
    input$target_plasmid,
    primer_paths$screening_path,
    "screening_offtarget_check.txt"
  )
  invisible(TRUE)
}

append_design_log <- function(log_path, stage, status, detail = "") {
  line <- paste(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stage,
    status,
    detail,
    sep = "\t"
  )
  cat(line, "\n", file = log_path, append = TRUE, sep = "")
}

run_design_stage <- function(stage, log_path, expression) {
  append_design_log(log_path, stage, "START")
  tryCatch(
    {
      value <- force(expression)
      append_design_log(log_path, stage, "OK")
      value
    },
    error = function(e) {
      reason <- conditionMessage(e)
      error_stage <- if (!is.null(e$stage)) e$stage else stage
      append_design_log(log_path, error_stage, "ERROR", reason)
      stop(structure(
        list(
          message = reason,
          call = NULL,
          stage = error_stage,
          parent = e
        ),
        class = c("design_stage_error", "error", "condition")
      ))
    }
  )
}

design_from_grna_pool <- function(
  input,
  feature,
  grnas,
  design_class,
  target_dir,
  log_path
) {
  attempted_ranges <- new.env(parent = emptyenv())
  attempts <- 0L
  last_failure_stage <- "homology_arms"
  result <- visit_grna_sets(
    grnas,
    input$parameters$n20_mn,
    input$parameters$n20_strands,
    function(candidates) {
      selected <- list(
        table = candidates,
        n20_range = range(c(candidates$n20_start, candidates$n20_end))
      )
      range_key <- paste(selected$n20_range, collapse = ":")
      if (exists(range_key, envir = attempted_ranges, inherits = FALSE)) {
        return(NULL)
      }
      assign(range_key, TRUE, envir = attempted_ranges)
      attempts <<- attempts + 1L
      detail <- sprintf(
        "set=%d;n20=%s;range=%s",
        attempts,
        paste(candidates$genomic_location, collapse = ","),
        range_key
      )
      append_design_log(log_path, "homology_arms", "TRY", detail)
      arms <- design_homology_arms(
        input,
        feature,
        selected,
        design_class,
        target_dir
      )
      if (is.null(arms)) {
        last_failure_stage <<- "homology_arms"
        append_design_log(
          log_path,
          "homology_arms",
          "REJECTED",
          sprintf("set=%d", attempts)
        )
        return(NULL)
      }
      append_design_log(
        log_path,
        "design_outputs",
        "TRY",
        sprintf("set=%d", attempts)
      )
      output_attempt <- tryCatch(
        list(
          ok = TRUE,
          value = write_design_outputs(
            input,
            feature,
            selected,
            arms,
            design_class,
            target_dir
          )
        ),
        error = function(e) list(ok = FALSE, error = e)
      )
      if (!output_attempt$ok) {
        last_failure_stage <<- "design_outputs"
        append_design_log(
          log_path,
          "design_outputs",
          "REJECTED",
          sprintf(
            "set=%d;reason=%s",
            attempts,
            conditionMessage(output_attempt$error)
          )
        )
        return(NULL)
      }
      append_design_log(
        log_path,
        "design_outputs",
        "OK",
        sprintf("set=%d", attempts)
      )
      list(
        selected = selected,
        arms = arms,
        primer_paths = output_attempt$value,
        attempts = attempts
      )
    }
  )
  if (is.null(result)) {
    reason <- sprintf(
      paste(
        "Не удалось завершить дизайн ни для одного допустимого набора N20",
        "(проверено уникальных диапазонов: %d)"
      ),
      attempts
    )
    stop(structure(
      list(
        message = reason,
        call = NULL,
        stage = last_failure_stage
      ),
      class = c("grna_sets_exhausted", "error", "condition")
    ))
  }
  result
}

write_target_error <- function(
  path,
  gene_name,
  design_class,
  stage,
  reason
) {
  writeLines(
    c(
      paste("timestamp", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "\t"),
      paste("gene", gene_name, sep = "\t"),
      paste("class", design_class, sep = "\t"),
      paste("stage", stage, sep = "\t"),
      paste("reason", reason, sep = "\t")
    ),
    path
  )
}

design_target <- function(input, genome_name, gene_name, design_class) {
  safe_name <- tolower(gsub("[^A-Za-z0-9_.-]", "_", gene_name))
  layout <- output_layout(input$output_dir)
  target_dir <- file.path(
    layout$tech_report,
    paste0(safe_name, "_results")
  )
  wet_lab_dir <- file.path(
    layout$wet_lab,
    paste0(safe_name, "_results")
  )
  unlink(file.path(
    wet_lab_dir,
    c(
      "final_sequences.fasta",
      "final_sequences.txt",
      "wet_lab_report.txt"
    )
  ))
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(target_dir, "design.log")
  error_path <- file.path(target_dir, "error.txt")
  writeLines(
    paste(
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      "target",
      "START",
      sprintf("gene=%s;class=%s", gene_name, design_class),
      sep = "\t"
    ),
    log_path
  )
  if (file.exists(error_path)) {
    unlink(error_path)
  }
  message(sprintf("[%s] %s", design_class, gene_name))
  tryCatch(
    {
      feature <- run_design_stage(
        "feature_lookup",
        log_path,
        feature_record(input, gene_name)
      )
      run_design_stage(
        "chopchop",
        log_path,
        run_chopchop(
          input,
          genome_name,
          feature,
          design_class,
          target_dir
        )
      )
      grnas <- run_design_stage(
        "n20_filter",
        log_path,
        filter_grnas(
          file.path(target_dir, "n20_table.tsv"),
          feature,
          design_class,
          input$parameters$n20_offtarget
        )
      )
      grnas <- run_design_stage(
        "n20_pool",
        log_path,
        prepare_grna_pool(
          grnas,
          input$parameters$n20_mn,
          input$parameters$n20_strands
        )
      )
      design <- run_design_stage(
        "homology_arms",
        log_path,
        design_from_grna_pool(
          input,
          feature,
          grnas,
          design_class,
          target_dir,
          log_path
        )
      )
      run_design_stage(
        "n20_output",
        log_path,
        write_selected_grnas(design$selected, feature, target_dir)
      )
      run_design_stage(
        "virtual_pcr",
        log_path,
        run_virtual_pcr(input, target_dir, design$primer_paths)
      )
      run_design_stage(
        "wet_lab_output",
        log_path,
        write_wet_lab_outputs(
          wet_lab_dir,
          feature,
          design_class,
          design$primer_paths$wet_lab$sequences,
          design$primer_paths$wet_lab$sequence_purposes,
          design$primer_paths$wet_lab$primer_metrics,
          design$primer_paths$wet_lab$screening_product_sizes
        )
      )
      append_design_log(log_path, "target", "OK")
      data.frame(
        gene = gene_name,
        class = design_class,
        output_dir = target_dir,
        status = "ok",
        stage = NA_character_,
        reason = NA_character_,
        wet_lab_dir = wet_lab_dir
      )
    },
    error = function(e) {
      stage <- if (inherits(e, "design_stage_error")) {
        e$stage
      } else {
        "unknown"
      }
      reason <- conditionMessage(e)
      write_target_error(
        error_path,
        gene_name,
        design_class,
        stage,
        reason
      )
      append_design_log(log_path, "target", "ERROR", paste(stage, reason))
      stop(structure(
        list(
          message = reason,
          call = NULL,
          gene = gene_name,
          design_class = design_class,
          stage = stage,
          target_dir = target_dir
        ),
        class = c("target_design_error", "error", "condition")
      ))
    }
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cli <- parse_designer_args(args)
  input <- make_design_input(cli)
  dir.create(input$output_dir, recursive = TRUE, showWarnings = FALSE)
  layout <- output_layout(input$output_dir)
  dir.create(layout$wet_lab, recursive = TRUE, showWarnings = FALSE)
  dir.create(layout$tech_report, recursive = TRUE, showWarnings = FALSE)
  targets <- bind_rows(
    data.frame(gene = cli$cds, class = rep("cds", length(cli$cds))),
    data.frame(gene = cli$ncrna, class = rep("ncrna", length(cli$ncrna)))
  )
  write_run_parameters(
    input,
    targets,
    file.path(layout$tech_report, "run_parameters.tsv")
  )
  genome_assets <- prepare_chopchop_assets(input)
  chopchop_config <- configure_chopchop(input, genome_assets)
  if (!file.copy(
    chopchop_config,
    file.path(layout$tech_report, "chopchop_config.json"),
    overwrite = TRUE
  )) {
    stop("Не удалось сохранить конфигурацию CHOPCHOP в TechReport", call. = FALSE)
  }
  results <- lapply(seq_len(nrow(targets)), function(i) {
    tryCatch(
      design_target(
        input,
        genome_assets$name,
        targets$gene[[i]],
        targets$class[[i]]
      ),
      error = function(e) {
        stage <- if (inherits(e, "target_design_error")) {
          e$stage
        } else {
          "unknown"
        }
        output_dir <- if (inherits(e, "target_design_error")) {
          e$target_dir
        } else {
          NA_character_
        }
        message(sprintf(
          "[ERROR] gene=%s class=%s stage=%s: %s",
          targets$gene[[i]],
          targets$class[[i]],
          stage,
          conditionMessage(e)
        ))
        data.frame(
          gene = targets$gene[[i]],
          class = targets$class[[i]],
          output_dir = output_dir,
          status = "error",
          stage = stage,
          reason = conditionMessage(e),
          wet_lab_dir = NA_character_
        )
      }
    )
  })
  write_tsv(
    bind_rows(results),
    file.path(layout$tech_report, "design_summary.tsv")
  )
}

if (sys.nframe() == 0) {
  main()
}
