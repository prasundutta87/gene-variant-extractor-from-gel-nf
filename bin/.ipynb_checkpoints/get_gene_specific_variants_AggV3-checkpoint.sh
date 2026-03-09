#!/bin/bash

ROW_INDEX="$1"
GENES_BED="$2"
BIALLELIC_SHARDS="$3"
ANNO_SHARDS="$4"
SITEQC_SHARDS="$5"

LINE=`head -${ROW_INDEX} ${GENES_BED}|tail -1`

printf "$LINE" >temp_${ROW_INDEX}.bed

chrom=`cat temp_${ROW_INDEX}.bed|cut -f1`
start=`cat temp_${ROW_INDEX}.bed|cut -f2`
end=`cat temp_${ROW_INDEX}.bed|cut -f3`
gene=`cat temp_${ROW_INDEX}.bed|cut -f4`

##get genotype information
bedtools intersect -wo -a temp_${ROW_INDEX}.bed -b ${BIALLELIC_SHARDS} >VCF_genotypes_${ROW_INDEX}.bed

##get functional annotation information
bedtools intersect -wo -a temp_${ROW_INDEX}.bed -b ${ANNO_SHARDS}>VCF_annotations_${ROW_INDEX}.bed

##get site QC information
bedtools intersect -wo -a temp_${ROW_INDEX}.bed -b ${SITEQC_SHARDS} >VCF_siteQC_${ROW_INDEX}.bed

#Add variant file header
printf "sample\tvariant_id\tFILTER\tgenotype\twhole_cohort_AC\twhole_cohort_AN\twhole_cohort_AF\n" > "$gene"_genotypes.tsv

#Add siteQC file header
printf "variant_id\tAC_Hom\tAC_Het\tAC_Hemi\tmedianDP\tmedianGQ\tABratio\tmissingness_rate\n" > "$gene"_siteQC.tsv

#Add annotation file header
first_vcf_line=$(cat VCF_annotations_${ROW_INDEX}.bed |head -n 1| cut -f13)
printf "variant_id\t$(bcftools +split-vep -l "$first_vcf_line" -l | cut -f 2 | paste -s -d '\t')\n" > "$gene"_annotation.tsv

##Convert 0-based BED start coordinate to 1-based for bcftools
region_start=$((start+1))

#Loop through multiple lines if a longer gene is present across two VCF chunks
cat VCF_genotypes_${ROW_INDEX}.bed| while read line;do
	variant_vcf=$(echo "$line" | cut -f13)
	bcftools view -Ou -r "$chrom":"$region_start"-"$end" -f 'PASS,.' --threads 4 "$variant_vcf" |  \
	bcftools filter -i 'AF <= 0.01' --threads 4 -Ou | \
	bcftools query -i 'GT="alt"' -f '[%SAMPLE\t%CHROM\_%POS\_%REF\_%ALT\t%FILTER\t%GT\t%AC\t%AN\t%AF\n]' >>"$gene"_genotypes.tsv
done

cat VCF_siteQC_${ROW_INDEX}.bed| while read line;do
	siteQC_vcf=$(echo "$line" | cut -f13)
	bcftools view -Ou -r "$chrom":"$region_start"-"$end" --threads 4 "$siteQC_vcf" |  \
	bcftools query -f '%CHROM\_%POS\_%REF\_%ALT\t%INFO/AC_Hom\t%AC_Het\t%AC_Hemi\t%MEDIAN_DP\t%MEDIAN_GQ\t%AB_RATIO\t%MISSINGNESS_RATE\n' >>"$gene"_siteQC.tsv
done

cat VCF_annotations_${ROW_INDEX}.bed| while read line;do
    annotation_vcf=$(echo "$line" | cut -f13)
	bcftools +split-vep -r "$chrom":"$region_start"-"$end" "$annotation_vcf" -d \
    -i "SYMBOL == \"$gene\" && (MANE_SELECT != \".\" || CANONICAL == \"YES\")" \
	-a CSQ -A tab \
	-f '%CHROM\_%POS\_%REF\_%ALT\t%CSQ\n' >>"$gene"_annotation.tsv
done

rm temp_${ROW_INDEX}.bed VCF_genotypes_${ROW_INDEX}.bed VCF_annotations_${ROW_INDEX}.bed VCF_siteQC_${ROW_INDEX}.bed

#Combine variants with their annotations
Rscript combine_variant_using_duckplyr.R $PWD "$gene"
echo "all steps successfully over for $gene"


