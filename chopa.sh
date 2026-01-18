#!/usr/bin/bash
eval "$(conda shell.bash hook)"
conda activate chopchop

# $1 - genome
# $2 - tsv
# $3 - gene name
# $4 - pTarget file (in .gb format)
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


TARGET_REGION=$(Rscript get_cords.R "$1" "$2" "$3")

# Проверяем, что TARGET_REGION не пуст
if [ -z "$TARGET_REGION" ]; then
  echo -e "Ошибка: Rscript не вернул координаты." >&2
  exit 1
fi


echo -e "Координаты таргета: $TARGET_REGION"

mkdir -p "$output_dir"

# Запуск chopchop
python valenlab-chopchop-a56388468523/chopchop.py \
  -Target "$TARGET_REGION" \
  -G "$genome_name" \
  -M NGG \
  -T 1 \
  -g 20 \
  -m 10 \
  --padSize 0 \
  --scoringMethod DOENCH_2016 \
  -o "$output_dir" > "$output_dir"/n20_table.tsv

mkdir -p "$output_dir"/offtargets_n20
mv "$output_dir"/*offtargets "$output_dir"/offtargets_n20/

echo -e "Извлекаю последовательности плеч гомологии..."
Rscript get_gRNA_place.R "$output_dir" "$1" "$4"

if [ $? -ne 0 ]; then
  echo -e "Ошибка при извлечении последовательностей плеч гомологии" >&2
  exit 1
fi


echo -e "\\n !!Проверка на оффтаргет для праймеров!! \\n"
java -jar virtualPCR/dist/virtualPCR.jar "$output_dir"/pcr_config.conf

rm "$output_dir"/sequence.fa "$output_dir"/gene_file.fa
mkdir -p "$output_dir"/additional_files
mv "$output_dir"/bowtie.err "$output_dir"/additional_files/
mv "$output_dir"/twoBitToFa.err "$output_dir"/additional_files/
mv "$output_dir"/right_arm_report.txt "$output_dir"/additional_files/
mv "$output_dir"/left_arm_report.txt "$output_dir"/additional_files/
mv "$output_dir"/pcr_config.conf "$output_dir"/additional_files/
mv "$output_dir"/primer3_settings.txt "$output_dir"/additional_files/
mv "$output_dir"/offtarget_check.txt "$output_dir"/additional_files/

mkdir -p "$output_dir"/sequences
mv "$output_dir"/offtarget.fa "$output_dir"/sequences/
mv "$output_dir"/homology_arms.fa "$output_dir"/sequences/

echo -e " \\nРезультаты сохранены в директории" "$output_dir"