#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(read.gb))

# 1 - result directory
# 2 - genome
# 2 - pTarget

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Требуется 3 аргумента: result_directory, genome", call. = FALSE)
}

result_dir <- args[1]
genome_path <- args[2]
pTarget_path <- args[3]


# Загружаем функцию callPrimer3 из GitHub
#suppressMessages(devtools::source_url("https://gist.githubusercontent.com/IdoBar/5e78ae7a5cc7277a04b126ce6f595d6e/raw/45c60662f3479f41765bce839835c4988a7e5b36/callPrimer3.R"))
source("callPrimer3.R")


# Загружаем геном
genome <- readDNAStringSet(genome_path, nrec=1)[[1]]
genome_length <- length(genome)

# Читаем таблицу с результатами ChopChop
gRNA_table_path <- file.path(result_dir, 'n20_table.tsv')
gRNA_table <- read_tsv(gRNA_table_path, show_col_types = FALSE) %>%
  janitor::clean_names() %>% 
  mutate(genomic_location = as.numeric(gsub('contig_\\d:', '', genomic_location)))

# Проверяем, что таблица не пустая
if (nrow(gRNA_table) == 0) {
  stop("Таблица с результатами ChopChop пустая", call. = FALSE)
}

# Выбираем лучшие gRNA, чтобы они были не на одной стрэнде
plus_strand <- gRNA_table %>% filter(strand == '+') %>% slice_head(n = 2)
minus_strand <- gRNA_table %>% filter(strand == '-') %>% slice_head(n = 2)
best_gRNA <- bind_rows(plus_strand, minus_strand) %>% 
  arrange(desc(efficiency)) %>%
  slice_head(n = 3)

if (nrow(best_gRNA) < 3) {
  stop("Недостаточно gRNA для определения границ вырезания (менее 2)", call. = FALSE)
}

# Определяем границы вырезания
range_borders <- 30
cut_region_range <- range(best_gRNA$genomic_location) + c(-range_borders, range_borders)
cut_region_range[1] <- max(1, cut_region_range[1])
cut_region_range[2] <- min(genome_length, cut_region_range[2])

cat(sprintf("Область вырезания: %d-%d\n", cut_region_range[1], cut_region_range[2]))

lha_length <- 300
rha_length <- 400
# Извлекаем последовательности плеч гомологии
left_arm_start <- max(1, cut_region_range[1] - (lha_length*1.5))
left_arm_end <- cut_region_range[1] - 1
right_arm_start <- cut_region_range[2] + 1
right_arm_end <- min(genome_length, cut_region_range[2] + (rha_length*1.5))

# Безопасное извлечение последовательностей
left_arm <- genome[left_arm_start:left_arm_end]
right_arm <- genome[right_arm_start:right_arm_end]

# Проверяем, что последовательности не пустые
if (length(left_arm) == 0) {
  stop(sprintf("Левое плечо пустое! Границы: %d-%d", left_arm_start, left_arm_end), call. = FALSE)
}
if (length(right_arm) == 0) {
  stop(sprintf("Правое плечо пустое! Границы: %d-%d", right_arm_start, right_arm_end), call. = FALSE)
}

# Создаем DNAStringSet для плеч гомологии
homology_arms <- DNAStringSet(list(
  left_arm_seq = left_arm,
  right_arm_seq = right_arm
))
names(homology_arms) <- c(
  paste0("left_arm_", left_arm_start, "_", left_arm_end),
  paste0("right_arm_", right_arm_start, "_", right_arm_end)
)

# Сохраняем плечи гомологии в FASTA
homology_arms_path <- file.path(result_dir, 'homology_arms.fa')
writeXStringSet(homology_arms, homology_arms_path)
homology_arms_path <- normalizePath(homology_arms_path)
#cat(sprintf("Последовательности плеч гомологии сохранены в: %s\n", homology_arms_path))



# === ВЫЗОВ PRIMER3 ДЛЯ ГЕНЕРАЦИИ ПРАЙМЕРОВ ===
cat("\n=== ГЕНЕРАЦИЯ ПРАЙМЕРОВ С ПОМОЩЬЮ PRIMER3 ===\n")

# Определяем пути к файлам и настройкам
primer3_path <- "primer3/src/primer3_core"
thermo_params_path <- "primer3/src/primer3_config/"

min_prod_size <- min(length(right_arm), length(left_arm))
max_prod_size <- max(length(right_arm), length(left_arm))
product_range <- c(min_prod_size, max_prod_size*1.5)

# Создаем файл настроек для primer3
settings_path <- file.path(result_dir, "primer3_settings.txt")
writeLines(c(
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
  "PRIMER_MIN_TM=50.0",
  "PRIMER_OPT_TM=60.0",
  "PRIMER_MAX_TM=65.0",
  "PRIMER_PAIR_MAX_DIFF_TM=8.0",
  "PRIMER_MIN_GC=40.0",
  "PRIMER_OPT_GC_PERCENT=50.0",
  "PRIMER_MAX_GC=60.0",
  "PRIMER_MAX_SELF_ANY=12.0",
  "PRIMER_MAX_SELF_END=8.0",
  "PRIMER_PAIR_MAX_COMPL_ANY=12.0",
  "PRIMER_PAIR_MAX_COMPL_END=8.0",
  "PRIMER_MAX_HAIRPIN_TH=47.0",
  "PRIMER_MAX_POLY_X=5",
  "PRIMER_MAX_NS_ACCEPTED=0",
  "PRIMER_SALT_CORRECTIONS=1",
  "PRIMER_DNA_CONC=50.0",
  "PRIMER_WT_SIZE_LT=0.5",
  "PRIMER_WT_SIZE_GT=0.5",
  "PRIMER_WT_TM_LT=0.5",
  "PRIMER_WT_TM_GT=0.5",
  "PRIMER_PAIR_WT_DIFF_TM=0.2",
  #  paste0("PRIMER_MISPRIMING_LIBRARY=", offtarget_path),
  #  "PRIMER_MAX_LIBRARY_MISPRIMING=15.0",
  #  "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=30.0",
  paste0("PRIMER_PRODUCT_SIZE_RANGE=", product_range[1], "-", product_range[2]),
  "PRIMER_EXPLAIN_FLAG=1",
  "PRIMER_FIRST_BASE_INDEX=1",
  "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0",
  "PRIMER_MAX_END_STABILITY=9.0",
  "PRIMER_PAIR_MAX_LIB_SIM=12.0",
  "PRIMER_LIBERAL_BASE=1",
  "="
), settings_path)


# Проверяем существование primer3_core
if (!file.exists(primer3_path)) {
  stop(sprintf("primer3_core не найден по пути: %s", primer3_path), call. = FALSE)
}

# Вызываем функцию callPrimer3 для левого плеча
cat("\nГенерирую праймеры для левого плеча...\n")
l_l <-c(lha_length, length(left_arm))
left_primers <- tryCatch({
  callPrimer3(
    seq = as.character(left_arm),
    size_range = paste0(l_l[1], '-', l_l[2]),
    Tm = c(50.0, 60.0, 65.0),
    Tm_diff = 10.0,
    name = names(homology_arms)[1],
    primer_num = 5,
    primer3 = primer3_path,
    thermo.param = thermo_params_path,
    report = file.path(result_dir, 'left_arm_report.txt'),
    settings = settings_path
  )
}, error = function(e) {
  warning(sprintf("Ошибка при генерации праймеров для левого плеча: %s", e$message))
  return(NA)
})
l_r <-c(rha_length, length(right_arm))
# Вызываем функцию callPrimer3 для правого плеча
cat("\nГенерирую праймеры для правого плеча...\n")
right_primers <- tryCatch({
  callPrimer3(
    seq = as.character(right_arm),
    size_range = paste0(l_r[1], '-', l_r[2]),
    Tm = c(57.0, 60.0, 63.0),
    Tm_diff = 5.0,
    name = names(homology_arms)[2],
    primer_num = 5,
    primer3 = primer3_path,
    thermo.param = thermo_params_path,
    report = file.path(result_dir, 'right_arm_report.txt'),
    settings = settings_path
  )
}, error = function(e) {
  warning(sprintf("Ошибка при генерации праймеров для правого плеча: %s", e$message))
  return(NA)
})

left_primers <- left_primers %>% mutate(genome_start=(PRIMER_LEFT_pos+left_arm_start), genome_end=(PRIMER_RIGHT_pos+left_arm_start))
right_primers <- right_primers %>% mutate(genome_start=(PRIMER_LEFT_pos+right_arm_start), genome_end=(PRIMER_RIGHT_pos+right_arm_start))
ha_primers <- bind_rows(left_primers, right_primers)
best_ha_primers <- ha_primers %>% filter(PRIMER_NUM == 1)
offtarget_range <- best_ha_primers %>% select(genome_start, genome_end) %>% range()

# Создаем последовательность для проверки офф-таргетов
offtarget_path <- file.path(result_dir, 'offtarget.fa')
offtg1 <- genome[1:offtarget_range[1]]
offtg2 <- genome[offtarget_range[2]:length(genome)]

DNAStringSet(list(offtarget1 = offtg1, offtarget2 = offtg2)) %>% writeXStringSet(offtarget_path)

#Скрининг генома праймеры
screening_range <- offtarget_range+c(-200, 200)
screening_seq <- genome[screening_range[1]:screening_range[2]]
genome_primers <- tryCatch({
  callPrimer3(
    seq = as.character(screening_seq),
    size_range = paste0(length(screening_seq)-100, '-', length(screening_seq)),
    Tm = c(50.0, 60.0, 65.0),
    Tm_diff = 5.0,
    name = paste0('genome_screening:', screening_range[1], '-', screening_range[2]),
    primer_num = 5,
    primer3 = primer3_path,
    thermo.param = thermo_params_path,
    report = file.path(result_dir, 'genome_screening_report.txt'),
    settings = settings_path
  )
}, error = function(e) {
  warning(sprintf("Ошибка при генерации праймеров для скрининга: %s", e$message))
  return(NA)
}) %>% 
  mutate(genome_start=(PRIMER_LEFT_pos+screening_range[1]), genome_end=(PRIMER_RIGHT_pos+screening_range[1]))
genome_screening_f <- genome_primers[1,3] %>% DNAString()
genome_screening_r <- genome_primers[1,4] %>% DNAString()
genome_screening_primers <- DNAStringSet(list(genome_screening_f,
                                         genome_screening_r))
names(genome_screening_primers) <- paste0('genome_screening_', c('F', 'R'))
genome_screening_tm_f <- genome_primers[1,5] %>% round(1)
genome_screening_tm_r <- genome_primers[1,6] %>% round(1)
genome_screening_product_size_before <- genome_primers[1,26]
genome_screening_start <- genome_primers[1,28]
genome_screening_old_end <- genome_primers[1,29]


#Результаты
gRNAs_forward_primers <- DNAStringSet(paste0('ACGACTAGT', 
                                             substring(best_gRNA$target_sequence, 1, 20), 
                                             'GTTTTAGAGCTAGAAATAGCAAGTTaaaataaggct'))
gRNA_revers_primer <- DNAStringSet(c(gRNA_revers = 'AGTTGACGCTAAAAAAAGCACCGACTCGGTGCC'))
names(gRNAs_forward_primers) <- paste0('gRNA_forward', 1:3)
left_homology_arm_primer <- DNAStringSet(c(paste0('AGCGTCAACT',left_primers$PRIMER_LEFT_SEQUENCE[1]),
                                           paste0('CTTGCGGGCAGTCAT',left_primers$PRIMER_RIGHT_SEQUENCE[1])))
names(left_homology_arm_primer) <- paste0('left_homology_arm_primer_', c('f', 'r'))
right_homology_arm_primer <- DNAStringSet(c(paste0('ATGACTGCCCGCAAG',right_primers$PRIMER_LEFT_SEQUENCE[1]), 
                                            paste0('ACGCTGCAG',right_primers$PRIMER_RIGHT_SEQUENCE[1])))
names(right_homology_arm_primer) <- paste0('right_homology_arm_primer_', c('f', 'r'))

c(gRNAs_forward_primers, 
  gRNA_revers_primer, 
  left_homology_arm_primer, 
  right_homology_arm_primer, 
  genome_screening_primers) %>%
  writeXStringSet(file.path(result_dir,'all_primers.fasta'))


virtual_pcr_config <- file.path(result_dir, "pcr_config.conf")
writeLines(c(
  paste0("targets_path=",offtarget_path),
  paste0("output_path=", file.path(result_dir, 'offtarget_check.txt')),
  paste0("primers_path=",file.path(result_dir, 'all_primers.fasta')),
  "type=primer",
  "ShowPCRProducts=true"
), virtual_pcr_config)


#pTarget res
pTargert <- read.gb(pTarget_path)$Exported

pTargert_seq <- DNAStringSet(c(pTarget_origin = pTargert$ORIGIN))[[1]]
ins_start <- regexpr('ACTAGT',pTargert_seq, ignore.case=T)[1]
ins_end <- regexpr('CTGCAG',pTargert_seq, ignore.case=T)[1]

pcr_products_gRNA <- sapply(gRNAs_forward_primers, function(x){
  c(x, DNAString('agtccgttatcaacttgaaaaagtggcaccgagtcggtgctttttttAGCGTCAACT'))
})
names(pcr_products_gRNA) <- paste0('gRNA_F', 1:3, '_pcr_product')

starts <- best_ha_primers$PRIMER_LEFT_pos + best_ha_primers$PRIMER_LEFT_len
ends <- best_ha_primers$PRIMER_RIGHT_pos - best_ha_primers$PRIMER_RIGHT_len

pcr_product_lha <- c(left_homology_arm_primer$left_homology_arm_primer_f, 
                     left_arm[starts[1]:ends[1]], 
                     (left_homology_arm_primer$left_homology_arm_primer_r %>% complement() %>% reverse()))
pcr_product_rha <- c(right_homology_arm_primer$right_homology_arm_primer_f, 
                     right_arm[starts[2]:ends[2]], 
                     (right_homology_arm_primer$right_homology_arm_primer_r %>% complement() %>% reverse()))

concatinate_pcr_result <- sapply(pcr_products_gRNA, function(x){
  c(x, pcr_product_lha[-(1:10)], pcr_product_rha[-(1:15)])
})
names(concatinate_pcr_result) <- paste0('concatinate_pcr_product_', 1:3)

PCR_product_sequnces <- DNAStringSet(c(pcr_products_gRNA, 
               'pcr_product_lha' = pcr_product_lha,
               'pcr_product_lha' = pcr_product_rha,
               concatinate_pcr_result))
writeXStringSet(PCR_product_sequnces, file.path(result_dir, 'PCR_product_sequnces.fasta'))

pTarget_sequences <- sapply(PCR_product_sequnces[6:8], function(x){
  c(pTargert_seq[1:ins_start],
    subseq(x, start = 4, end = (length(x) - 3)),
    pTargert_seq[ins_end:length(pTargert_seq)])
}) %>% DNAStringSet()
names(pTarget_sequences) <- paste0('pTarget_', 1:3)

writeXStringSet(pTarget_sequences, file.path(result_dir, 'pTargets.fasta'))

edited_genome <- DNAStringSet(c(edited_genome = c(offtg1,
                                                  left_arm[starts[1]:ends[1]],
                                                  right_homology_arm_primer$right_homology_arm_primer_f,
                                                  right_arm[starts[2]:ends[2]],
                                                  offtg2)))
writeXStringSet(edited_genome, file.path(result_dir, 'edited_genome.fasta'))

pTarget_start <- regexpr('GCATCTGTGCGGTATTTCAC', pTargert_seq)[1]
old_end <- (DNAString('TGCTTATGGAGCTGCACATG') %>% 
              complement() %>% 
              reverse() %>% 
              regexpr(pTargert_seq))
new_end <- (DNAString('TGCTTATGGAGCTGCACATG') %>% 
              complement() %>% 
              reverse() %>% 
              regexpr(pTarget_sequences[[1]]))
primer_screning_before <- (old_end[1]+attr(old_end, "match.length")) - pTarget_start
primer_screning_after <- (new_end[1]+attr(new_end, "match.length")) - pTarget_start

genome_screening_new_end <- genome_screening_r %>% 
  complement() %>%
  reverse() %>%
  regexpr(edited_genome)
genome_screening_product_size_after <- (genome_screening_new_end[1]+attr(genome_screening_new_end, "match.length")) - genome_screening_start




report_path <- file.path(result_dir, 'report.txt')
writeLines(c(
  paste0('Температура отжига прямого праймера для левого плеча гомологии=', round(best_ha_primers[1,5], 1)),
  paste0('Температура отжига обратного праймера для левого плеча гомологии=', round(best_ha_primers[1,6], 1)),
  paste0('Температура отжига прямого праймера для правого плеча гомологии=', round(best_ha_primers[2,5], 1)),
  paste0('Температура отжига обратного праймера для правого плеча гомологии=', round(best_ha_primers[2,6], 1)),
  '',
  paste0('Длина продукта скрининга плазмиды до модификации=', primer_screning_before),
  paste0('Длина продукта скрининга плазмиды после модификации=', primer_screning_after),
  paste0('Длина продукта скрининга генома до модификации=', genome_screening_product_size_before),
  paste0('Длина продукта скрининга генома после модификации=', genome_screening_product_size_after),
  '',
  paste0('Температура отжига прямого праймера для скрининга генома=', genome_screening_tm_f),
  paste0('Температура отжига обратного праймера для скрининга генома=', genome_screening_tm_r)
), report_path)
