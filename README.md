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

## SoMoSeq single-nucleus RNA-seq analysis

-   R Markdown notebooks for the analysis of the SoMoSeq dataset (PIK3CA- and SLC35A2-related focal MCD cases and neurotypical autopsy controls).

-   These notebooks use the zUMIs derived gene expression matrix (intron+exon counts), from the generated FASTQ-to-count preprocessing workflow (`../scripts/` dir) and the sample- and cell-level metadata deposited on Zenodo.

-   Data: [Zenodo record](https://zenodo.org/records/22693814)

## Notebooks

| Notebook | Step | Writes (`work/`; deposited object/file name) |
|------------------|-----------------------|-------------------------------|
| `01_QC_preprocessing.Rmd` | Sample annotations, sample selection, nucleus QC thresholds, gene filters | ``` 01_qc_filtered.rds;``RDS/filtered_annotated_object.rds ``` (counts) |
| `02_normalization_clustering.Rmd` | Library-size log normalization, batch-blocked HVGs, HVG PCA, SNN/Leiden clustering | `02_clustered.rds`, ``` 02_sex_chromosome_flags.tsv;``normalized_expression/ ```, `hvg_pca` in the RDS |
| `03_celltype_annotation.Rmd` | Marker-contrast scoring and gated assignment | `03_annotated.rds`, ``` 03_annotation_per_cell.tsv;``metadata/final_cell_type_annotations.tsv.gz ``` |
| `04_annotation_validation.Rmd` | Coverage, marker enrichment, held-out panel, score–depth associations, SingleR agreement | ``` 04_annotation_validation/;``metadata/reference_annotations/ ``` |
| `05_pseudobulk_DE.Rmd` | Donor x batch x genotype x cell-type pseudobulk, limma-voom | ``` 05_differential_expression/;``results/pseudobulk/ ```, `results/differential_expression/` |
| `06_pathway_analysis.Rmd` | fgsea on moderated t ranks | ``` 06_pathways/;``results/pathways/ ``` |

Each notebook checks the final output against the deposited files. Notebooks 04–06 can also use the deposited annotated object as a starting point.

## Setup

1.  Download and extract the Zenodo folders. The notebooks expect the `SoMoSeq/` dir, which contains `RDS/`, `metadata/` and `results/`.
2.  Download the gene sets (below).
3.  Copy `config/example_paths.yml` to `config/paths.yml` and edit the paths. Relative paths are resolved from `config/`.
4.  Render the notebooks from `notebooks/`, for example `rmarkdown::render("01_QC_preprocessing.Rmd")`.

### Input data structure

-   `RDS/raw_unfiltered_object.rds`: a `SingleCellExperiment` of 44,716 genes x 1,116 nuclei, with raw `counts`.
    -   rowData: `symbol`, `gene_biotype`, `ensembl_id`.
    -   colData used by notebook 01: `cell_id`, `Sample_name`, `Seq_batch`, `Genotype`, `Case_control`, `Bulk_sc`, `Parent_cDNA_plate_ID`, `is_duplicate_group_member`, `spliced_counts`, `total_counts`.

### Gene sets

-   [MSigDB v2024.1 (human, gene symbols)](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/)

| File                            | Sets used |
|---------------------------------|-----------|
| `h.all.v2024.1.Hs.symbols.gmt`  | 50        |
| `c2.cp.v2024.1.Hs.symbols.gmt`  | 18        |
| `c5.hpo.v2024.1.Hs.symbols.gmt` | 24        |
| `c6.all.v2024.1.Hs.symbols.gmt` | 6         |

-   Epilepsy / cortical-malformation HPO sets (`epilepsy_gmt`; 13 sets).
    -   The HPO sets come from the [MSigDB gene-set pages](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/%3CNAME%3E):
        -   `HP_ABNORMAL_CEREBRAL_CORTEX_MORPHOLOGY`
        -   `HP_ABNORMAL_CEREBRAL_VASCULAR_MORPHOLOGY`
        -   `HP_CEREBRAL_CAVERNOUS_MALFORMATION`
        -   `HP_CORTICAL_DYSPLASIA`
        -   `HP_EEG_WITH_FOCAL_EPILEPTIFORM_DISCHARGES`
        -   `HP_EPILEPTIC_ENCEPHALOPATHY`
        -   `HP_FOCAL_CORTICAL_DYSPLASIA_TYPE_II`
        -   `HP_FOCAL_CORTICAL_DYSPLASIA`
        -   `HP_HEMIMEGALENCEPHALY`
        -   `HP_MEGALENCEPHALY`
        -   `HP_MULTIFOCAL_EPILEPTIFORM_DISCHARGES`.
-   The file also contains two sets from [Enrichr libraries](https://maayanlab.cloud/Enrichr/#libraries):
    -   "Proteins Involved in Epilepsy", from the `Elsevier_Pathway_Collection` library;
    -   "Cortical Dysplasia With Other Brain Malformations", from the `GeDiPNet_2023` library.

### SingleR references

-   This notebook uses the deposited SingleR predictions

-   To rerun SingleR, set `run_singler: true` in notebook 04 and point `config/paths.yml` to reference objects built from:

-   siletti:

    -   Siletti et al. (2023), "Transcriptomic diversity of cell types across the adult human brain", *Science*, <https://doi.org/10.1126/science.add7046>
    -   Cell type col = `cluster_annotation`.

-   developmental:

    -   Zhu et al. (2023), "Multi-omic profiling of the developing human cerebral cortex at the single-cell level", *Science Advances*, <https://doi.org/10.1126/sciadv.adg3754>
    -   10x multiome (RNA + ATAC) from 12 donors, two at each of six ages: 18–19 GW, 23–24 GW, 0 y, 4–6 y, 14 y and 20–39 y
    -   Cell type col = `cluster_annotation`.

-   ss2:

    -   Allen Institute human MTG Smart-seq2 (Hodge et al. 2019, *Nature*)
    -   Cell type col = `subclass_label`.

## Software versions

-   R 4.5.2
-   scuttle 1.20.0
-   scran 1.38.0
-   scater 1.38.0
-   igraph 2.2.1
-   edgeR 4.8.2
-   limma 3.66.0
-   fgsea 1.39.4
-   SingleR 2.10.0
