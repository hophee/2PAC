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
    annotation <- read_tsv(path, skip = skip_lines, show_col_types = FALSE) %>%
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
      "Неподдерживаемый формат аннотации. Допустимы: bakta, gff",
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
    help = "Annotation format: bakta or gff"
  )
  parser <- add_argument(
    parser,
    "--chopchop-script",
    help = "Path to chopchop.py"
  )
  parser <- add_argument(parser, "--primer3", help = "Path to primer3_core")
  parser <- add_argument(
    parser,
    "--virtual-pcr-jar",
    help = "Path to virtualPCR.jar"
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
    targets[nzchar(targets)]
  }

  parsed <- parse_args(parser, args)
  annotation_format <- normalize_scalar(parsed$annotation_format)
  if (!length(annotation_format)) {
    annotation_format <- "gff"
  }
  annotation_format <- tolower(annotation_format)
  if (!annotation_format %in% c("bakta", "gff")) {
    stop(
      "Некорректный --annotation-format. Допустимы: bakta, gff",
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
    virtual_pcr_jar = normalize_scalar(parsed$virtual_pcr_jar)
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

split_gene_list <- function(values) {
  if (!length(values)) {
    return(character())
  }
  genes <- trimws(unlist(strsplit(values, ",", fixed = TRUE)))
  unique(genes[nzchar(genes)])
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

make_design_input <- function(cli) {
  annotation_format <- if (length(cli$annotation_format)) {
    cli$annotation_format[[1]]
  } else {
    "bakta"
  }
  # TODO: only a complete single-contig bacterial genome is supported for now.
  # Multi-contig FASTA/GFF input needs contig-aware sequence extraction.
  genome_set <- readDNAStringSet(cli$genome[[1]], nrec = 1)
  list(
    genome_path = cli$genome[[1]],
    genome = genome_set[[1]],
    genome_contig = names(genome_set)[[1]],
    annotation = read_genome_annotation(
      cli$genome_annotation[[1]],
      annotation_format
    ),
    target_plasmid = cli$target_plasmid[[1]],
    cas_plasmid = if (length(cli$cas_plasmid)) cli$cas_plasmid[[1]] else NULL,
    output_dir = cli$output_dir[[1]],
    tools = list(
      chopchop_script = if (length(cli$chopchop_script)) {
        cli$chopchop_script[[1]]
      } else {
        "chopchop/chopchop.py"
      },
      primer3 = if (length(cli$primer3)) {
        cli$primer3[[1]]
      } else {
        "primer3/src/primer3_core"
      },
      primer3_config = "primer3/src/primer3_config",
      virtual_pcr_jar = if (length(cli$virtual_pcr_jar)) {
        cli$virtual_pcr_jar[[1]]
      } else {
        "virtualPCR/dist/virtualPCR.jar"
      }
    )
  )
}

feature_record <- function(input, gene_name) {
  idx <- find_target_feature(input$annotation, gene_name)
  row <- input$annotation[idx, , drop = FALSE]
  contig <- if ("seqid" %in% names(input$annotation)) {
    as.character(row$seqid[[1]])
  } else {
    NA_character_
  }
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
  genome_dir <- dirname(normalizePath(input$genome_path))
  two_bit <- file.path(genome_dir, paste0(genome_name, ".2bit"))
  bowtie_prefix <- file.path(genome_dir, genome_name)
  if (!file.exists(two_bit)) {
    run_tool("faToTwoBit", c(input$genome_path, two_bit))
  }
  if (!file.exists(paste0(bowtie_prefix, ".1.ebwt"))) {
    run_tool("bowtie-build", c(input$genome_path, bowtie_prefix))
  }
  list(name = genome_name, directory = genome_dir)
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

select_grnas <- function(table_path, feature, design_class, target_dir) {
  grnas <- read_tsv(table_path, show_col_types = FALSE) %>%
    janitor::clean_names() %>%
    mutate(
      genomic_start = as.numeric(sub("^.*:", "", genomic_location)),
      genomic_end = genomic_start + 22
    ) %>%
    mutate(
      mid_closeness = abs(
        (feature$start + feature$end) %/%
          2 -
          (genomic_start + genomic_end) %/% 2
      ) /
        feature$length
    ) %>%
    filter(self_complementarity == 0, mm0 == 0)
  if (design_class == "cds") {
    grnas <- filter(grnas, mid_closeness <= 0.18)
  }
  candidates <- bind_rows(
    filter(grnas, strand == "+") %>% slice_head(n = 2),
    filter(grnas, strand == "-") %>% slice_head(n = 2)
  ) %>%
    arrange(mid_closeness) %>%
    slice_head(n = 3)
  if (nrow(candidates) < 3) {
    stop("Недостаточно подходящих N20 (требуется 3)", call. = FALSE)
  }
  write_tsv(candidates, file.path(target_dir, "selected_n20_table.tsv"))
  n20 <- DNAStringSet(substr(candidates$target_sequence, 1, 20))
  names(n20) <- paste0(feature$display_name, "_n20_", seq_along(n20))
  writeXStringSet(n20, file.path(target_dir, "selected_n20.fasta"))
  list(
    table = candidates,
    range = range(c(candidates$genomic_start, candidates$genomic_end)) +
      c(-40, 40)
  )
}

write_primer3_settings <- function(path, left_length, right_length) {
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

choose_pair <- function(primers, plus_strand, feature, design_class) {
  left <- primers$left
  right <- primers$right
  if (!nrow(left) || !nrow(right)) {
    return(NULL)
  }
  if (plus_strand) {
    deleted <- outer(right$genome_start, left$genome_end, "-") - 1
  } else {
    deleted <- outer(left$genome_start, right$genome_end, "-") - 1
  }
  eligible <- deleted < feature$length
  # TODO: add the legacy fallback that shortens the bridge when no in-frame
  # primer pair exists; the current CDS implementation requires %in% 3 == 0.
  if (design_class == "cds") {
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
  left_start <- if (feature$strand == "+") {
    max(1L, interval[[1]] - 350L)
  } else {
    interval[[2]] + 1L
  }
  left_end <- if (feature$strand == "+") {
    interval[[1]] - 1L
  } else {
    min(length(input$genome), interval[[2]] + 350L)
  }
  right_start <- if (feature$strand == "+") {
    interval[[2]] + 1L
  } else {
    max(1L, interval[[1]] - 450L)
  }
  right_end <- if (feature$strand == "+") {
    min(length(input$genome), interval[[2]] + 450L)
  } else {
    interval[[1]] - 1L
  }
  settings <- file.path(target_dir, "primer3_settings.txt")
  write_primer3_settings(
    settings,
    abs(left_end - left_start) + 1,
    abs(right_end - right_start) + 1
  )
  if (!file.exists(input$tools$primer3)) {
    stop("primer3_core не найден", call. = FALSE)
  }

  for (attempt in seq_len(1000)) {
    if (feature$strand == "+") {
      left_end <- interval[[1]] - 1L
      right_start <- interval[[2]] + 1L
      left_seq <- input$genome[left_start:left_end]
      right_seq <- input$genome[right_start:right_end]
    } else {
      right_end <- interval[[1]] - 1L
      left_start <- interval[[2]] + 1L
      left_seq <- complement(reverse(input$genome[left_start:left_end]))
      right_seq <- complement(reverse(input$genome[right_start:right_end]))
    }
    if (!length(left_seq) || !length(right_seq)) {
      stop("Плечо гомологии пусто", call. = FALSE)
    }
    writeXStringSet(
      DNAStringSet(c(left_arm = left_seq, right_arm = right_seq)),
      file.path(target_dir, "homology_arms_before_primer_search.fasta")
    )
    tm <- if (design_class == "cds") c(62.5, 63, 63.5) else c(60.5, 61, 62.5)
    left <- tryCatch(
      callPrimer3(
        as.character(left_seq),
        paste0("300-", length(left_seq)),
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
        paste0("400-", length(right_seq)),
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
        design_class
      )
      if (!is.null(pair)) {
        ticks <- sort(unlist(pair[, c("genome_start", "genome_end")]))
        if (ticks[[2]] >= feature$start && ticks[[3]] <= feature$end) {
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
    if (
      interval[[1]] > selected$range[[1]] ||
        interval[[2]] < selected$range[[2]] ||
        diff(interval) + 1 < feature$length * .3
    ) {
      break
    }
  }
  stop("Не удалось подобрать праймеры для плеч гомологии", call. = FALSE)
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
  frame_status <- "not_applicable"
  if (design_class == "cds") {
    bridge <- substr(bridge, 1, 15 - (gap %% 3))
    frame_status <- if (gap %% 3 == 0) {
      "in-frame"
    } else {
      sprintf("bridge shortened to %d nt", nchar(bridge))
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

  writeLines(
    c(
      paste("target", feature$query_name, sep = "\t"),
      paste("class", design_class, sep = "\t"),
      paste("deleted_nt", gap, sep = "\t"),
      paste("frame_status", frame_status, sep = "\t"),
      paste("left_primer_tm", round(pair$PRIMER_LEFT_TM[[1]], 1), sep = "\t"),
      paste("right_primer_tm", round(pair$PRIMER_RIGHT_TM[[2]], 1), sep = "\t"),
      paste(
        "screening_forward_tm",
        round(screening$PRIMER_LEFT_TM[[1]], 1),
        sep = "\t"
      ),
      paste(
        "screening_reverse_tm",
        round(screening$PRIMER_RIGHT_TM[[1]], 1),
        sep = "\t"
      )
    ),
    file.path(target_dir, "report.tsv")
  )
  list(
    all_primers_path = file.path(target_dir, "all_primers.fasta"),
    plain_path = plain_path,
    screening_path = screening_path
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

design_target <- function(input, genome_name, gene_name, design_class) {
  feature <- feature_record(input, gene_name)
  safe_name <- tolower(gsub("[^A-Za-z0-9_.-]", "_", gene_name))
  target_dir <- file.path(input$output_dir, paste0(safe_name, "_results"))
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("[%s] %s", design_class, gene_name))
  run_chopchop(input, genome_name, feature, design_class, target_dir)
  selected <- select_grnas(
    file.path(target_dir, "n20_table.tsv"),
    feature,
    design_class,
    target_dir
  )
  arms <- design_homology_arms(
    input,
    feature,
    selected,
    design_class,
    target_dir
  )
  primer_paths <- write_design_outputs(
    input,
    feature,
    selected,
    arms,
    design_class,
    target_dir
  )
  run_virtual_pcr(input, target_dir, primer_paths)
  data.frame(
    gene = gene_name,
    class = design_class,
    output_dir = target_dir,
    status = "ok"
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cli <- parse_designer_args(args)
  input <- make_design_input(cli)
  dir.create(input$output_dir, recursive = TRUE, showWarnings = FALSE)
  genome_assets <- prepare_chopchop_assets(input)
  configure_chopchop(input, genome_assets)
  cds <- split_gene_list(cli$cds)
  ncrna <- split_gene_list(cli$ncrna)
  targets <- bind_rows(
    data.frame(gene = cds, class = rep("cds", length(cds))),
    data.frame(gene = ncrna, class = rep("ncrna", length(ncrna)))
  )
  results <- lapply(seq_len(nrow(targets)), function(i) {
    tryCatch(
      design_target(
        input,
        genome_assets$name,
        targets$gene[[i]],
        targets$class[[i]]
      ),
      error = function(e) {
        message(sprintf(
          "[ERROR] %s (%s): %s",
          targets$gene[[i]],
          targets$class[[i]],
          conditionMessage(e)
        ))
        data.frame(
          gene = targets$gene[[i]],
          class = targets$class[[i]],
          output_dir = NA_character_,
          status = conditionMessage(e)
        )
      }
    )
  })
  write_tsv(
    bind_rows(results),
    file.path(input$output_dir, "design_summary.tsv")
  )
}

if (sys.nframe() == 0) {
  main()
}
