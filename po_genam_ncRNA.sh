#!/usr/bin/env bash
if [ -z "$1" ]; then
  echo "Usage: $0 <output_directory>"
  exit 1
fi

# Формируем имя файла лога на основе первого аргумента
LOG_FILE="log_${1}.txt"
rm "$LOG_FILE"
echo -e 'Лог пишется в файл ' "$(readlink -f $LOG_FILE) \n"

# Перенаправляем весь вывод в лог
exec > >(tee -a "$LOG_FILE") 2>&1

eval "$(conda shell.bash hook)"

mkdir -p "$1"
for eczv_id in $(cat ncRNAs.txt); do
    ./chopa_ncRNA.sh ./ZVL2_anot/ZvL2_Glu_clone_1_hn_re.fna ./ZVL2_anot/ZvL2_Glu_clone_1_hn_re.tsv "$eczv_id" plasmids.fasta "$1"/"$eczv_id"_results
done

cd "$1" || { echo "Failed to cd into $1"; exit 1; }
zip "$1"_results.zip *_results/*
mkdir -p collected

find . -type f -path './*_results/sequences/homology_arms.fasta' -print0 | \
while IFS= read -r -d '' f; do
    d=$(basename "$(dirname "$(dirname "$f")")")
    cp -n -- "$f" "collected/${d}__offtargets.fasta"
done

cat -- collected/*.fasta > all_ha_"$1".fasta