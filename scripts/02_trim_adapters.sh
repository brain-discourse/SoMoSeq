#!/bin/bash
# ==========================================================
# SoMoSeq preprocessing - step 2: adapter trimming (cutadapt)
# ==========================================================
#SBATCH -t 10-00:00:00
#SBATCH --mem=200g
#SBATCH -n 20

# usage: sbatch 02_trim_adapters.sh <RAW_DATA_DIR>
# <RAW_DATA_DIR> must contain reads_for_zUMIs.R1.fastq.gz and reads_for_zUMIs.R2.fastq.gz
# (produced by step 1, merging demultiplexed lane/sample fastqs).
# Trims SmartSeq-2/SoMoSeq TSO and Illumina adapter sequences from combined R1/R2 fastqs
# using cutadapt's linked-adapter/anchored adapter strategy (separate -g/-a per end), 
# quality-trims ends at Phred 15, and 
# drops reads shorter than 20bp after trimming. See paths.example.yaml for RAW_DATA_DIR and README.md for the full pipeline this script sits in.

module load cutadapt

RAW_DATA_DIR=$1
if [[ -z "$RAW_DATA_DIR" ]]; then
    echo "Usage: $0 <RAW_DATA_DIR>"
    exit 1
fi

cutadapt \
    -g AAGCAGTGGTATCAACGCAGAGTAC \
    -g CATGGAAGCAGTGGTATCAACGCAGAGTAC \
    -g ATGGGAAGCAGTGGTATCAACGCAGAGTAC \
    -a AAGCAGTGGTATCAACGCAGAGTACATGGA \
    -a AAGCAGTGGTATCAACGCAGAGTACATGGG \
    -a AAGCAGTGGTATCAACGCAGAGTACTTTTT \
    -G AAGCAGTGGTATCAACGCAGAGTAC \
    -G CATGGAAGCAGTGGTATCAACGCAGAGTAC \
    -G ATGGGAAGCAGTGGTATCAACGCAGAGTAC \
    -A AAGCAGTGGTATCAACGCAGAGTACATGGA \
    -A AAGCAGTGGTATCAACGCAGAGTACATGGG \
    -A AAGCAGTGGTATCAACGCAGAGTACTTTTT \
    --minimum-length=20 \
    --overlap=15 \
    -e 0.1 -q 15,15 \
    --trim-n --cores=20 \
    -o $RAW_DATA_DIR/reads_for_zUMIs_trimmed.re.R1.fastq.gz \
    -p $RAW_DATA_DIR/reads_for_zUMIs_trimmed.re.R2.fastq.gz \
    $RAW_DATA_DIR/reads_for_zUMIs.R1.fastq.gz $RAW_DATA_DIR/reads_for_zUMIs.R2.fastq.gz > $RAW_DATA_DIR/reads_for_zUMIs_trim.re.log
