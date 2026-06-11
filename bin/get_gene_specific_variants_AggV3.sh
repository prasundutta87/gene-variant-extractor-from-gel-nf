#!/bin/bash
# get_gene_specific_variants_AggV3.sh
# Extracts rare variants for a single gene from GEL AggV3 VCF shards.
# Pulls genotypes, siteQC metrics and VEP/GreenDB annotations from separate
# shard types, then merges them via R/duckplyr into a single per-gene TSV.
#
# Args (passed by Nextflow):
#   $1 genotype_vcfs   : file listing genotype shard paths (one per line)
#   $2 annotation_vcfs : file listing VEP/GreenDB annotation shard paths
#   $3 siteqc_vcfs     : file listing siteQC shard paths
#   $4 gene_bed_file   : single-row BED (chrom, start, end, gene_name)
#
# Outputs (to Nextflow working dir; staged by publishDir):
#   ${gene}_genotypes.tsv, ${gene}_siteQC.tsv, ${gene}_annotation.tsv, ${gene}.tsv
#
# Notes:
#   - chrX siteQC shards lack DP/GQ/ABratio/missingness — handled separately.
#   - While-read loops append across shards when a gene spans a boundary.
#   - No mkdir/mv — Nextflow publishDir handles output staging.

set -euo pipefail

genotype_vcfs="$1"
annotation_vcfs="$2"
siteqc_vcfs="$3"
gene_bed_file="$4"

chrom=$(cat "$gene_bed_file"|cut -f1)
start=$(cat "$gene_bed_file"|cut -f2)
end=$(cat "$gene_bed_file"|cut -f3)
gene=$(cat "$gene_bed_file"|cut -f4)

#Convert 0-based BED start coordinate to 1-based for bcftools
region_start=$((start+1))

echo "Processing: $gene ($chrom:$region_start-$end)"

#Add variant file header
printf "sample\tvariant_id\tFILTER\tgenotype\n" > "$gene"_genotypes.tsv

#Add siteQC file header
if [[ "$chrom" == "chrX" ]]; then
	printf "variant_id\tGEL_cohort_AC\tGEL_cohort_AN\tGEL_cohort_AF\tAC_Hom\tAC_Het\tAC_Hemi\tMEDIAN_DP_XX\tMEDIAN_DP_XY\tMEDIAN_GQ_XX\tMEDIAN_GQ_XY\tAB_RATIO_XX\tAB_RATIO_XY\tMISSINGNESS_RATE_XX\tMISSINGNESS_RATE_XY\n" > "$gene"_siteQC.tsv
else
	printf "variant_id\tGEL_cohort_AC\tGEL_cohort_AN\tGEL_cohort_AF\tAC_Hom\tAC_Het\tAC_Hemi\tmedianDP\tmedianGQ\tABratio\tmissingness_rate\n" > "$gene"_siteQC.tsv
fi

#Add annotation file header
first_vcf_line=$(cat "$annotation_vcfs" |head -n 1)
printf "variant_id\t$(bcftools +split-vep -l "$first_vcf_line" | cut -f 2 | paste -s -d '\t')\n" > "$gene"_annotation.tsv

#Loop through multiple lines if a longer gene is present across two VCF chunks
while read line;do
	variant_vcf="$line"
	bcftools view -Ou -r "$chrom":"$region_start"-"$end" -f 'PASS,.' --threads 4 "$variant_vcf" |  \
	bcftools query -i 'GT="alt"' -f '[%SAMPLE\t%CHROM\_%POS\_%REF\_%ALT\t%FILTER\t%GT\n]' >>"$gene"_genotypes.tsv
done < "$genotype_vcfs"

echo "All genotypes from all samples extracted"

if [[ "$chrom" == "chrX" ]]; then
	while read line;do
		siteQC_vcf="$line"
		bcftools view -Ou -r "$chrom":"$region_start"-"$end" --threads 4 "$siteQC_vcf" |  \
		bcftools query -f '%CHROM\_%POS\_%REF\_%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/AC_Hom\t%INFO/AC_Het\t%INFO/AC_Hemi\t%MEDIAN_DP_XX\t%MEDIAN_DP_XY\t%MEDIAN_GQ_XX\t%MEDIAN_GQ_XY\t%AB_RATIO_XX\t%AB_RATIO_XY\t%MISSINGNESS_RATE_XX\t%MISSINGNESS_RATE_XY\n' >>"$gene"_siteQC.tsv
	done < "$siteqc_vcfs"
else
	while read line;do
		siteQC_vcf="$line"
		bcftools view -Ou -r "$chrom":"$region_start"-"$end" --threads 4 "$siteQC_vcf" |  \
		bcftools query -f '%CHROM\_%POS\_%REF\_%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/AC_Hom\t%INFO/AC_Het\t%INFO/AC_Hemi\t%MEDIAN_DP\t%MEDIAN_GQ\t%AB_RATIO\t%MISSINGNESS_RATE\n' >>"$gene"_siteQC.tsv
	done < "$siteqc_vcfs"
fi

echo "All siteQC annotations extracted"

while read line;do
	annotation_vcf="$line"
	bcftools +split-vep -r "$chrom":"$region_start"-"$end" "$annotation_vcf" -d \
	-i "SYMBOL == \"$gene\" && (MANE_SELECT != \".\" || CANONICAL == \"YES\")" \
	-a CSQ -A tab \
	-f '%CHROM\_%POS\_%REF\_%ALT\t%CSQ\n' >>"$gene"_annotation.tsv
done < "$annotation_vcfs"

echo "All GreenDB and VEP annotations extracted"

#Combine sample genotypes, siteQC, GreenDB and VEP annotations
combine_variant_using_duckplyr.R "$PWD" "$gene"

echo "All steps successfully completed for $gene"