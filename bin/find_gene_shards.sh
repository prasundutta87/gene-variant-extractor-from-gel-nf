#!/bin/bash
# find_gene_shards.sh
#
# Finds overlapping AggV3 VCF shards for each gene in the input BED file.
# Called by Nextflow for each row in the genes BED file.
#
# NOTE: This script is not meant to be run directly.
#       Run the pipeline via Nextflow:
#           nextflow run main.nf --genes_bed my_genes.bed
#
#       Only --genes_bed needs to be specified - all other parameters
#       (shard manifests, reference files, output directory) are set
#       as defaults in nextflow.config.

ROW_INDEX="$1"
GENES_BED="$2"
BIALLELIC_GENOTYPE_SHARDS="$3"
ANNO_SHARDS="$4"
SITEQC_SHARDS="$5"

LINE=$(head -${ROW_INDEX} "${GENES_BED}" | tail -1)
echo "$LINE" > temp_${ROW_INDEX}.bed

gene=$(cut -f4 temp_${ROW_INDEX}.bed)

##get genotype information
bedtools intersect -wo -a temp_${ROW_INDEX}.bed -b "${BIALLELIC_GENOTYPE_SHARDS}" > VCF_genotypes_${ROW_INDEX}.bed

##get functional annotation information
bedtools intersect -wo -a temp_${ROW_INDEX}.bed -b "${ANNO_SHARDS}" > VCF_annotations_${ROW_INDEX}.bed

##get site QC information
bedtools intersect -wo -a temp_${ROW_INDEX}.bed -b "${SITEQC_SHARDS}" > VCF_siteQC_${ROW_INDEX}.bed

# Extract only the VCF paths from column-13
cat VCF_genotypes_${ROW_INDEX}.bed | cut -f13 > ${gene}_biallelic_genotype_shards.txt
cat VCF_annotations_${ROW_INDEX}.bed | cut -f13 > ${gene}_anno_shards.txt
cat VCF_siteQC_${ROW_INDEX}.bed | cut -f13 > ${gene}_siteqc_shards.txt

# Extract only the VCF index paths from column-14
cat VCF_genotypes_${ROW_INDEX}.bed | cut -f14 > ${gene}_biallelic_genotype_shards_idx.txt
cat VCF_annotations_${ROW_INDEX}.bed | cut -f14 > ${gene}_anno_shards_idx.txt
cat VCF_siteQC_${ROW_INDEX}.bed | cut -f14 > ${gene}_siteqc_shards_idx.txt

rm VCF_genotypes_${ROW_INDEX}.bed VCF_annotations_${ROW_INDEX}.bed VCF_siteQC_${ROW_INDEX}.bed
