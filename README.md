This repository contains a DSL2 Nextflow pipeline for gene-level variant extraction and aggregation.

The workflow:

  1) Iterates over regions defined in a BED file

  2) Extracts gene-specific rare variants (Agg V3 whole cohort AF <= 0.01) from biallelic shard/subshard datasets

  3) Integrates functional annotation and site-level QC information from their respective shards/subshard datasets

  4) Produces gene-level summary TSV outputs

  5) Executes reproducibly inside a Docker container

  6) The pipeline is designed for cloud execution and parallel processing, with one process per gene region.

All computational tools (bcftools, bedtools, samtools, R, duckplyr) are encapsulated within a Docker image to ensure reproducibility and portability.

