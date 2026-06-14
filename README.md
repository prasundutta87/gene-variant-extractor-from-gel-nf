# gene-variant-extractor-from-gel-nf

A DSL2 Nextflow pipeline for extracting, annotating, and classifying rare variants at the gene level from Genomics England (GEL) AggV3 germline short variant calls.

GEL's AggV3 dataset brings together germline genomes from 100kGP, NHS GMS, and COVID-19 participants, generated using Illumina DRAGEN's Iterative GVCF Genotyper (DRAGEN 3.7.8). Because of the sheer scale of the data, variants are split across genomic segments called "shards", with separate shard types for genotypes, site-level QC metrics, and functional annotation. This pipeline was built to make working with those shards tractable at the gene level. Give it a BED file of genes you care about and it handles the rest.

## How to run

```bash
nextflow run main.nf --genes_bed my_genes.bed
```

`--genes_bed` is the only parameter you need to provide. Everything else - shard manifests, reference files, phenotype data, output directory, is already set in `nextflow.config` and rarely needs changing outside of a data release update.

The BED file should be tab-separated with four columns and no header: chromosome, start (0-based), end, gene name. For example:

```
chr2    197386005    197408815    SF3B1
chrX    53225828     53321350     IQSEC2
```

Coordinates should be GRCh38. The gene name must match the SYMBOL field in the VEP annotation.

## What it does

The pipeline runs three processes per gene:

**FIND_SHARDS** : takes the input BED file and uses `bedtools intersect` to find which AggV3 VCF shards overlap each gene. Genes that span shard boundaries are handled automatically.

**RUN_GENE** : the main extraction step. Uses `bcftools` to pull genotypes, site-level QC metrics, and VEP/GreenDB functional annotations from their respective shards, then merges everything into a single per-gene TSV using `duckplyr`. Variants with a cohort AF > 1% are dropped at this stage. chrX genes are handled separately because GEL's siteQC reports sex-stratified metrics (XX/XY) for the X chromosome rather than the single-value metrics used for autosomes.

**ANNOTATE_GENE** : joins the per-gene variant table with participant metadata and RD/GMS phenotype data, then runs variant classification. This includes AF compatibility checks against gnomAD and the GEL cohort, ClinVar status, and per-variant call quality flags (depth, genotype quality, missingness). Outputs a fully annotated table per gene in both `tsv.gz` and `parquet` format.

## Outputs

Results land in `results/<GENE>/` for each gene:

- `<gene>_for_review.tsv.gz` : the annotated variant table, one row per participant per variant
- `<gene>_for_review.parquet` : the same table in parquet format, useful for faster querying across many genes

## Docker environment

All the tools the pipeline needs are bundled into a single Docker image (`prasundutta87/gene-variant-extractor-from-gel-docker-image:2.0.0`), so there's nothing to install manually. The image is based on `rocker/tidyverse:4.5.2` and includes:

- `bcftools 1.21` (compiled from source with S3/libcurl support for direct access to GEL's S3 buckets)
- `samtools`, `bedtools`, `tabix`
- R packages: `tidyverse`, `data.table`, `duckdb`, `duckplyr`, `arrow`, `gt`, `ggpubr`

The container is specified in `nextflow.config` and is pulled automatically when the pipeline runs. You just need Docker enabled on your executor.

## A note on the environment

This pipeline is designed to run on Lifebit CloudOS with access to GEL's AWS S3 data buckets. The shard paths, phenotype file locations, and other reference paths in `nextflow.config` all point to GEL-internal resources and won't be accessible outside the GEL research environment.

## Requirements

- Nextflow ≥ 25.10
- Docker
- Access to GEL AggV3 data via Lifebit CloudOS
