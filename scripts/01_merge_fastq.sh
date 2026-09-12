#!/bin/bash
# ==========================================================
# SoMoSeq preprocessing - step 1: combine demultiplexed fastqs
# ==========================================================
#SBATCH -t 10-00:00:00
#SBATCH --mem=200g
#SBATCH -N 1
#SBATCH -n 10

# usage: sbatch 01_merge_fastq.sh <raw_fastq_dir>
# <raw_fastq_dir> contains the per-sample demultiplexed fastqs for one lane/sample pool.
# Produces one concatenated R1 and R2 fastq (reads_for_zUMIs.R1/R2.fastq.gz), 
# an index fastq of cell-specific barcode sequences, and 
# a text file mapping each cell to its barcode using (merge_demultiplexed_fastq.R)[https://github.com/sdparekh/zUMIs/blob/main/misc/merge_demultiplexed_fastq.R]

module load r
echo 'R_LIBS_USER="${R_LIBS_USER}"' > $HOME/.Renviron
module load pigz

RAW_FASTQ_DIR=$1
if [[ -z "$RAW_FASTQ_DIR" ]]; then
    echo "Usage: $0 <raw_fastq_dir>"
    exit 1
fi

Rscript ${ZUMIS_INSTALL_DIR}/misc/merge_demultiplexed_fastq.R --dir $RAW_FASTQ_DIR
