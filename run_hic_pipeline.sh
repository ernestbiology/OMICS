#!/usr/bin/env bash
set -euo pipefail

############################################
### КОНФИГ ###
############################################

BASE_URL="https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics"

GENOME_SHORT="T2T_human"
ENZYME="DpnII"
THREADS=4

REF_DIR="data/reference"
REF_FASTA="$REF_DIR/${GENOME_SHORT}.fna"
CHROM_SIZES="$REF_DIR/chrom.sizes"
RESTRICTION_SITES="$REF_DIR/restriction_sites_${ENZYME}.txt"

JUICER_DIR="$(pwd)/tools/juicer"

declare -A R1_REMOTE=(
  [MoPh7]="Copy%20of%20MoPh7_S85_L001_R1_001.fastq.gz"
  [MoPh11]="Copy%20of%20MoPh11_S86_L001_R1_001.fastq.gz"
  [MoPh14]="Copy%20of%20MoPh14_S87_L001_R1_001.fastq.gz"
  [MoPh15]="Copy%20of%20MoPh15_S88_L001_R1_001.fastq.gz"
)
declare -A R2_REMOTE=(
  [MoPh7]="Copy%20of%20MoPh7_S85_L001_R2_001.fastq.gz"
  [MoPh11]="Copy%20of%20MoPh11_S86_L001_R2_001.fastq.gz"
  [MoPh14]="Copy%20of%20MoPh14_S87_L001_R2_001.fastq.gz"
  [MoPh15]="Copy%20of%20MoPh15_S88_L001_R2_001.fastq.gz"
)

SAMPLES=("MoPh7" "MoPh11" "MoPh14" "MoPh15")

############################################
### Директории ###
############################################

mkdir -p data/raw data/trimmed "$REF_DIR" data/juicer
mkdir -p results/fastqc_raw results/cutadapt results/hic
mkdir -p scripts tools

############################################
### БЛОК 1. Подготовка референса (один раз) ###
############################################

if [[ ! -f "$REF_FASTA" ]]; then
    echo "=== Скачивание и подготовка референса ==="

    wget -O "$REF_DIR/T2T_human.fna.gz" \
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"

    gzip -dkf "$REF_DIR/T2T_human.fna.gz"

    python3 scripts/rename_chroms_t2t.py

    bwa index "$REF_FASTA"
    samtools faidx "$REF_FASTA"
    cut -f1,2 "${REF_FASTA}.fai" > "$CHROM_SIZES"
else
    echo "=== Референс уже подготовлен, пропускаю ==="
fi

if [[ ! -d tools/juicer ]]; then
    echo "=== Установка Juicer ==="
    git clone \
      --branch juicer_course_version \
      --single-branch \
      https://github.com/dpanc2/OMICS_course_spring_2026.git \
      tools/juicer
fi

if [[ ! -f "$RESTRICTION_SITES" ]]; then
    echo "=== Генерация файла сайтов рестрикции ($ENZYME) ==="
    python3 tools/juicer/misc/generate_site_positions.py \
      "$ENZYME" "$GENOME_SHORT" "$REF_FASTA"
    mv "${GENOME_SHORT}_${ENZYME}.txt" "$RESTRICTION_SITES"
fi

############################################
### БЛОК 2. Цикл по образцам ###
############################################

for SAMPLE in "${SAMPLES[@]}"; do
    echo "=========================================="
    echo "=== Образец: $SAMPLE ==="
    echo "=========================================="

    RAW_R1="data/raw/${SAMPLE}_R1.fastq.gz"
    RAW_R2="data/raw/${SAMPLE}_R2.fastq.gz"

    if [[ ! -s "$RAW_R1" ]]; then
        wget --no-check-certificate -O "$RAW_R1" \
          "${BASE_URL}/${R1_REMOTE[$SAMPLE]}"
    fi
    if [[ ! -s "$RAW_R2" ]]; then
        wget --no-check-certificate -O "$RAW_R2" \
          "${BASE_URL}/${R2_REMOTE[$SAMPLE]}"
    fi

    fastqc "$RAW_R1" "$RAW_R2" -o results/fastqc_raw

    TRIMMED_R1="data/trimmed/${SAMPLE}_R1.trimmed.fastq.gz"
    TRIMMED_R2="data/trimmed/${SAMPLE}_R2.trimmed.fastq.gz"
    CUTADAPT_LOG="results/cutadapt/${SAMPLE}.cutadapt.log"

    cutadapt \
      -q 20 \
      -m 70 \
      -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
      -o "$TRIMMED_R1" \
      -p "$TRIMMED_R2" \
      "$RAW_R1" "$RAW_R2" \
      > "$CUTADAPT_LOG" 2>&1

    echo "cutadapt лог: $CUTADAPT_LOG"

    JUICER_SAMPLE_DIR="$(pwd)/data/juicer/${SAMPLE}"
    mkdir -p "$JUICER_SAMPLE_DIR/fastq"

    ln -sf "$(pwd)/$TRIMMED_R1" "$JUICER_SAMPLE_DIR/fastq/${SAMPLE}_R1.fastq.gz"
    ln -sf "$(pwd)/$TRIMMED_R2" "$JUICER_SAMPLE_DIR/fastq/${SAMPLE}_R2.fastq.gz"

    bash "$JUICER_DIR/scripts/juicer.sh" \
      -D "$JUICER_DIR" \
      -d "$JUICER_SAMPLE_DIR" \
      -g "$GENOME_SHORT" \
      -z "$(pwd)/$REF_FASTA" \
      -p "$(pwd)/$CHROM_SIZES" \
      -y "$(pwd)/$RESTRICTION_SITES" \
      -s "$ENZYME" \
      -t "$THREADS" \
      2>&1 | tee "$JUICER_SAMPLE_DIR/juicer.log"

    HIC_SRC="$JUICER_SAMPLE_DIR/aligned/inter_30.hic"
    HIC_DST="results/hic/${SAMPLE}.inter_30.hic"

    if [[ -f "$HIC_SRC" ]]; then
        cp "$HIC_SRC" "$HIC_DST"
        echo "OK: $HIC_DST"
    else
        echo "ОШИБКА: $HIC_SRC не найден, см. $JUICER_SAMPLE_DIR/juicer.log" >&2
        exit 1
    fi

done

echo "=========================================="
echo "Готовые .hic файлы:"
ls -la results/hic/*.hic
