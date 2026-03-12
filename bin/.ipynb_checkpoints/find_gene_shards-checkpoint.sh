#!/bin/bash

ROW_INDEX="$1"
GENES_BED="$2"
BIALLELIC_SHARDS="$3"
ANNO_SHARDS="$4"
SITEQC_SHARDS="$5"

LINE=$(head -${ROW_INDEX} "${GENES_BED}" | tail -1)
echo "$LINE" > temp_${ROW_INDEX}.bed

gene=$(cut -f4 temp_${ROW_INDEX}.bed)

##get genotype information
bedtools intersect -wo -a temp_${ROW_INDEX}.bed -b "${BIALLELIC_SHARDS}" > VCF_genotypes_${ROW_INDEX}.bed

##get functional annotation information
bedtools intersect -wo -a temp_${ROW_INDEX}.bed -b "${ANNO_SHARDS}" > VCF_annotations_${ROW_INDEX}.bed

##get site QC information
bedtools intersect -wo -a temp_${ROW_INDEX}.bed -b "${SITEQC_SHARDS}" > VCF_siteQC_${ROW_INDEX}.bed

# Extract only the VCF paths from column-13
cat VCF_genotypes_${ROW_INDEX}.bed | cut -f13 > ${gene}_biallelic_shards.txt
cat VCF_annotations_${ROW_INDEX}.bed | cute -f13 > ${gene}_annotation_shards.txt
cat VCF_siteQC_${ROW_INDEX}.bed |cut -f13 > ${gene}_siteqc_shards.txt

rm VCF_genotypes_${ROW_INDEX}.bed VCF_annotations_${ROW_INDEX}.bed VCF_siteQC_${ROW_INDEX}.bed
