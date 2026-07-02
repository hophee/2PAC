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
  janitor::clean_names() %>%
  filter(!is.na(gene))

if (grepl("ECZV_", args[3]) | grepl("ZvL2_Glu_", args[3])) {
  g_idx <- which(genome_table$locus_tag == args[3])[1]
} else {
  g_idx <- which(genome_table$gene == args[3])[1]
}
gene_start <- genome_table$start[g_idx]
gene_end <- genome_table$stop[g_idx]

gene_length <- (gene_end - gene_start + 1)

# ranges_get <- function(x){
#   if (x<=500){
#     return(x)
#   } else if(x>500 & x< 1500){
#     return(500+((x-500)/10))
#   } else{
#     return(x/3)
#   }
# }
#
# plot(x=1:3000, y=sapply(1:3000, ranges_get))

if (gene_length <= 500) {
  gap <- round(gene_length / 10, -1) * 2
  cut_start <- gene_start - gap
  cut_end <- gene_end + gap

  nearest_genes <- genome_table %>%
    slice(g_idx - 1, g_idx + 1) %>%
    select(start, stop) %>%
    as.matrix() %>%
    sort()
  nearest_genes <- nearest_genes[2:3]
  cut_start <- max(cut_start, nearest_genes[1])
  cut_end <- min(cut_end, nearest_genes[2])
} else if (gene_length > 500 & gene_length < 1500) {
  side <- 250
  mid <- gene_start + gene_length %/% 2
  cut_start <- mid - side
  cut_end <- mid + side
} else {
  cut_start <- floor((gene_end - gene_start) / 3) + gene_start
  cut_end <- gene_end - floor((gene_end - gene_start) / 3)
}
cat(paste0(genome_name, ":", cut_start, "-", cut_end))
