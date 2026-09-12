#!/bin/bash
# ==========================================================
# SoMoSeq preprocessing - step 0: pre-process QC + first-pass adapter removal (fastp)
# ==========================================================
#SBATCH -t 5-00:00:00
#SBATCH -n 4
#SBATCH --mem=80g

# usage: sbatch 00_preprocess_qc.sh <R1.fastq.gz> <R2.fastq.gz> <OUT_DIR>
# Runs on raw demultiplexed per-sample fastqs (provided by the sequencing core)
# Deduplicates, auto-detects and trims generic Illumina adapters, and separates out reads without a mate. 
# this is the first, generic adapter-removal step
# a second, SmartSeq-2/SoMoSeq TSO-specific pass happens later in scripts/02_trim_adapters.sh,

module load fastp
module load fastqc
module load multiqc

R1_FILE=$1
R2_FILE=$2
OUT_DIR=$3
if [[ -z "$R1_FILE" || -z "$R2_FILE" || -z "$OUT_DIR" ]]; then
    echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> <OUT_DIR>"
    exit 1
fi

SAMPLE_NAME=$(basename "$R1_FILE" _R1_001.fastq.gz)

fastp \
    --in1 "$R1_FILE" --in2 "$R2_FILE" \
    --out1 "$OUT_DIR/${SAMPLE_NAME}_dedup_R1_001.fastq.gz" \
    --out2 "$OUT_DIR/${SAMPLE_NAME}_dedup_R2_001.fastq.gz" \
    --dedup --overrepresentation_analysis -c --detect_adapter_for_pe --thread 4 \
    --failed_out "$OUT_DIR/${SAMPLE_NAME}_failed.fastq.gz" \
    --unpaired1 "$OUT_DIR/${SAMPLE_NAME}_unpaired_R1_001.fastq.gz" \
    --unpaired2 "$OUT_DIR/${SAMPLE_NAME}_unpaired_R2_001.fastq.gz"

# Run once per lane/sample pool after all samples are processed:
# fastqc $OUT_DIR/*.fastq.gz -o $OUT_DIR
# multiqc $OUT_DIR -o $OUT_DIR/multiqc_report
