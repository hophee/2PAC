#!/usr/bin/bash
#source: https://github.com/hophee/2PAC

eval "$(conda shell.bash hook)"
conda activate chopchop

# $1 - genome
# $2 - tsv
# $3 - gene name
# $4 - pTarget+pCas fasta
# $5 - output dir

output_dir=${5,,}
genome_name=$(basename "$1" | cut -d. -f1)

# 2bit
if [ ! -f "tbit/${genome_name}.2bit" ]; then
  mkdir -p tbit
  faToTwoBit "$1" "tbit/${genome_name}.2bit"
fi

# Проверка bowtie индексов
if [ ! -f "bwt_idx/${genome_name}.1.ebwt" ]; then
  mkdir -p bwt_idx
  bowtie-build "$1" "bwt_idx/${genome_name}"
fi


TARGET_REGION=$(Rscript get_cords2.R "$1" "$2" "$3")

# Проверяем, что TARGET_REGION не пуст
if [ -z "$TARGET_REGION" ]; then
  echo -e "Ошибка: Rscript не вернул координаты." >&2
  exit 1
fi


echo -e "Координаты таргета: $TARGET_REGION"

mkdir -p "$output_dir"

# Запуск chopchop
echo 'Запуск chopchop'
python chopchop/chopchop.py \
  -Target "$TARGET_REGION" \
  -G "$genome_name" \
  -M NGG \
  -T 1 \
  -g 20 \
  -m 3 \
  --padSize 0 \
  -O 50 \
  --scoringMethod DOENCH_2016 \
  -o "$output_dir" > "$output_dir"/n20_table.tsv

mkdir -p "$output_dir"/offtargets_n20
mv "$output_dir"/*offtargets "$output_dir"/offtargets_n20/

echo -e "Извлекаю последовательности плеч гомологии..."
Rscript get_gRNA_place2.R "$output_dir" "$1" "$4" "$2" "$3"

if [ $? -ne 0 ]; then
  echo -e "Ошибка при извлечении последовательностей плеч гомологии" >&2
  exit 1
fi


echo -e "\\n Проверка на оффтаргет для праймеров"
java -jar virtualPCR/dist/virtualPCR.jar "$output_dir"/pcr_config.conf &> /dev/null
java -jar virtualPCR/dist/virtualPCR.jar "$output_dir"/pcr_config2.conf &> /dev/null
java -jar virtualPCR/dist/virtualPCR.jar "$output_dir"/pcr_config3.conf &> /dev/null

rm "$output_dir"/sequence.fa "$output_dir"/gene_file.fa
mkdir -p "$output_dir"/additional_files
mv "$output_dir"/bowtie.err "$output_dir"/additional_files/
mv "$output_dir"/twoBitToFa.err "$output_dir"/additional_files/
mv "$output_dir"/right_arm_report.txt "$output_dir"/additional_files/
mv "$output_dir"/left_arm_report.txt "$output_dir"/additional_files/
mv "$output_dir"/pcr_config.conf "$output_dir"/additional_files/
mv "$output_dir"/pcr_config2.conf "$output_dir"/additional_files/
mv "$output_dir"/primer3_settings.txt "$output_dir"/additional_files/
mv "$output_dir"/genome_screening_report.txt "$output_dir"/additional_files/
mv "$output_dir"/primer3_table.tsv "$output_dir"/additional_files/

mkdir -p "$output_dir"/sequences
mv "$output_dir"/offtarget.fa "$output_dir"/sequences/
mv "$output_dir"/homology_arms.fasta "$output_dir"/sequences/
mv "$output_dir"/homology_arms_before_primer_search.fa "$output_dir"/sequences/

cp "$output_dir"/all_primers.fasta "$output_dir"/all_primers.txt

echo -e " \\n \\n!!Результаты сохранены в директории" "$(readlink -f $output_dir)!! \n"
