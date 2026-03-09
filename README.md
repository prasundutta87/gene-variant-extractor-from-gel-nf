This repository contains a DSL2 Nextflow pipeline for gene-level variant extraction and aggregation.

The workflow:

  Iterates over regions defined in a BED file

  Extracts gene-specific rare variants (Agg V3 whole cohort AF <= 0.01) from biallelic shard/subshard datasets

  Integrates functional annotation and site-level QC information from their respective shards/subshard datasets

  Produces gene-level summary TSV outputs

  Executes reproducibly inside a Docker container

  The pipeline is designed for cloud execution and parallel processing, with one process per gene region.

All computational tools (bcftools, bedtools, samtools, R, duckplyr) are encapsulated within a Docker image to ensure reproducibility and portability.
