suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# 1 - genome
# 2 - tsv
# 3 locus tag/gene name

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Требуется 3 аргумента: геном, tsv-таблица, название гена/локус тэг",
        call. = FALSE
    )
}
genome_name <- readDNAStringSet(args[1])@ranges@NAMES[1]


lines <- readLines(args[2])
hash_lines <- grep("^#", lines)
skip_lines <- if (length(hash_lines) > 0) max(hash_lines) - 1 else 0
genome_table <- read_tsv(args[2], skip = skip_lines, show_col_types = FALSE) %>%
    janitor::clean_names()

if (grepl("ECZV_", args[3])) {
    g_idx <- which(genome_table$locus_tag == args[3])[1]
} else {
    g_idx <- which(genome_table$gene == args[3])[1]
}
start <- genome_table$start[g_idx]
end <- genome_table$stop[g_idx]


if ((end - start + 1) > 500) {
    mid <- floor((end - start) / 2) + start
    side <- (((end - start) / 10 + 500) / 2) %>% floor()
    start <- mid - side + 1
    end <- mid + side
}
cat(paste0(genome_name, ":", start, "-", end))
