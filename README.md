# SoMoSeq

A method for genotype informed single nuclei RNA sequencing of mosaic brain tissue

## Preprocessing workflow

Raw FASTQs are processed through five steps to produce gene-expression count matrices,
run once per sample / lane/sample pool (see `sample_sheet.tsv`):

0. **`scripts/00_preprocess_qc.sh`**
   - first-pass QC and generic Illumina adapter removal on demultiplexed per-sample fastqs
   - uses fastp (dedup, `--detect_adapter_for_pe`, filter reads without their mate), followed by FastQC/MultiQC.
   - Run per sample, before samples are merged into a lane/sample pool.
1. **`scripts/01_merge_fastq.sh`**
   - combines per-sample demultiplexed FASTQs for a lane/sample pool into one concatenated R1/R2 fastq via zUMIs `merge_demultiplexed_fastq.R`,
   - also generates a cell-barcode index fastq and
   - a barcode-to-cell map
2. **`scripts/02_trim_adapters.sh`**
   - a second, targeted adapter-trimming step
   - run on the combined R1/R2 FASTQs with cutadapt
   - trims SmartSeq-2/SoMoSeq TSO adapter sequences (min length 20bp, Phred 15 quality trim, 10% error rate, 15bp minimum
     overlap)
   - not the same as step_00 which is strictly generic adapter trimming
3. **`scripts/03_filter_index.pl`**
   - filters the barcode index fastq file (from `scripts/01_merge_fastq.sh` ) to only keep pairs that survived targeted adapter-trimming, preserving read order
4. **`scripts/04_run_zumis.sh`**
   - runs [zUMIs](https://github.com/sdparekh/zUMIs)(version 2.9.7e) which performs the following against a per-lane/sample config (see `configs/zUMIs_template.run.yaml`):
     4.1. _BC filtering_ - uses Phred 20 - cells with fewer than 100 reads are dropped
     4.2. _STAR alignment_ - aligns against a GRCh38/GENCODE v32 STAR index
     4.3. _feature counting_ - includes both exonic and intronic reads - primary hits retained only for multimappers

Each step's script applies identical logic across all lanes/samples by directory/config path (see `sample_sheet.tsv` for the per-lane/sample values)

## Configuration

- Copy `configs/paths.example.yaml`, fill in the values for your environment, and substitute them into `configs/zUMIs_template.run.yaml` and the `scripts/*.sh` calls before running.
- See the header comments in each script/config for expected inputs and outputs.
- write to @brain-discourse for exact per-lane/sample parameter values actually used for this project (implemented on UNC's longleaf cluster via slurm workflow manager)

## Dependencies

- fastp v0.23.4
- FastQC v0.12.1
- MultiQC v1.28
- cutadapt (v 4.9 with Python 3.9.6),
- zUMIs (v 2.9.7e with zUMIs internal miniconda env), includes STAR (v 2.7.3a), Rsubread (v 1.32.4), samtools (v 1.9), Perl (`IO::Zlib`), R (v 3.6.0)

## Notes

### GRCh38_GENCODE reference used for STAR index

- Obtained Primary sequence from GRCh38 from GENCODE (11/22/19)
  `wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz`
  `wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz`

- Created two files with mapping of Gene/Transript ID's to gene names
  `awk -F '\t|\"' '($3 ~ /transcript/){print $12,$16}' gencode.v32.primary_assembly.annotation.gtf > Ens_transcript_ID_to_Gene_Name`
  `awk -F '\t|\"' '($3 ~ /gene/){print $10,$14}' gencode.v32.primary_assembly.annotation.gtf > Ens_Gene_ID_to_Gene_Name`
