This repository contains a DSL2 Nextflow pipeline for gene-level variant extraction and aggregation from Genomics England (GEL)'s Aggregated Variant Calls (AggV3). AggV3 brings together short variants in germline genomes from 100kGP, NHS GMS and Covid-19 participants. AggV3 was prepared with by Illumina DRAGEN's Iterative GVCF Genotyper using genomes aligned using the DRAGEN 3.7.8 pipeline. Due to the size of the data, there are multiple VCFs, each representing a segment of the genome, known as "shards" and "subshards".

The workflow:

  1) Iterates over regions defined in a BED file

  2) Extracts gene-specific rare variants (Agg V3 whole cohort AF <= 0.01) from biallelic shard/subshard datasets

  3) Integrates functional annotation and site-level QC information from their respective shards/subshard datasets

  4) Produces gene-level summary TSV outputs

  5) Executes reproducibly inside a Docker container

  6) The pipeline is designed for cloud execution and parallel processing, with one process per gene region.

All computational tools (bcftools, bedtools, samtools, R, duckplyr) are encapsulated within a Docker image to ensure reproducibility and portability.

This pipeline is intended for systematic gene-level variant characterization and supports integration of annotation and QC metrics into structured tabular outputs suitable for downstream analysis.



