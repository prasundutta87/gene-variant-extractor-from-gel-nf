#!/bin/bash

genotype_vcfs="$1"
annotation_vcfs="$2"
siteqc_vcfs="$3"
gene_bed_file="$4"

chrom=`cat "$gene_bed_file"|cut -f1`
start=`cat "$gene_bed_file"|cut -f2`
end=`cat "$gene_bed_file"|cut -f3`
gene=`cat "$gene_bed_file"|cut -f4`

##Convert 0-based BED start coordinate to 1-based for bcftools
region_start=$((start+1))

#Add variant file header
printf "sample\tvariant_id\tFILTER\tgenotype\n" > "$gene"_genotypes.tsv

#Add siteQC file header
printf "variant_id\tGEL_cohort_AC\tGEL_cohort_AN\tGEL_cohort_AF\tAC_Hom\tAC_Het\tAC_Hemi\tmedianDP\tmedianGQ\tABratio\tmissingness_rate\n" > "$gene"_siteQC.tsv

#Add annotation file header
first_vcf_line=$(cat "$annotation_vcfs" |head -n 1)
printf "variant_id\t$(bcftools +split-vep -l "$first_vcf_line" -l | cut -f 2 | paste -s -d '\t')\n" > "$gene"_annotation.tsv

#Loop through multiple lines if a longer gene is present across two VCF chunks
while read line;do
	variant_vcf="$line"
	bcftools view -Ou -r "$chrom":"$region_start"-"$end" -f 'PASS,.' --threads 4 "$variant_vcf" |  \
	bcftools filter -i 'AF <= 0.01' --threads 4 -Ou | \
	bcftools query -i 'GT="alt"' -f '[%SAMPLE\t%CHROM\_%POS\_%REF\_%ALT\t%FILTER\t%GT\t%AC\t%AN\t%AF\n]' >>"$gene"_genotypes.tsv
done < "$genotype_vcfs"

while read line;do
	siteQC_vcf="$line"
	bcftools view -Ou -r "$chrom":"$region_start"-"$end" --threads 4 "$siteQC_vcf" |  \
	bcftools query -f '%CHROM\_%POS\_%REF\_%ALT\t%AC\t%AN\t%AF\t%AC_Hom\t%AC_Het\t%AC_Hemi\t%MEDIAN_DP\t%MEDIAN_GQ\t%AB_RATIO\t%MISSINGNESS_RATE\n' >>"$gene"_siteQC.tsv
done < "$siteqc_vcfs"

while read line;do
    annotation_vcf="$line"
	bcftools +split-vep -r "$chrom":"$region_start"-"$end" "$annotation_vcf" -d \
    -i "SYMBOL == \"$gene\" && (MANE_SELECT != \".\" || CANONICAL == \"YES\")" \
	-a CSQ -A tab \
	-f '%CHROM\_%POS\_%REF\_%ALT\t%CSQ\n' >>"$gene"_annotation.tsv
done < "$annotation_vcfs"

#Combine variants with their annotations
Rscript combine_variant_using_duckplyr.R "$PWD" "$gene"
echo "all steps successfully over for $gene"


