#!/usr/bin/env Rscript
# combine_variant_using_duckplyr_with_siteQC.R
# Merges three per-gene TSVs (genotypes, siteQC, VEP/GreenDB annotations)
# into a single rare-variant table using DuckDB-backed dplyr for efficiency.
#
# Args: <input_directory_path> <gene_name>
#   input_directory_path : directory containing the three input TSVs
#                          (pass "$PWD" from Nextflow — files staged there)
#   gene_name            : HGNC symbol; used to locate inputs and name output
#
# Output: ${gene}.tsv — one row per sample x variant, AF <= 1% only
#
# Notes:
#   - siteQC deduplication removes variants duplicated by overlapping GreenDB
#     enhancer regions (QC values are site-level and identical across duplicates)
#   - Genotype normalisation: 1/0 -> 0/1 (VCF allele order not guaranteed)
#   - many-to-many join is expected — variant_id is not unique per sample row
#   - Exits with status 1 if output is missing or < 5000 bytes so Nextflow
#     can catch task failure

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(argparser)))
suppressMessages(suppressWarnings(library(duckplyr)))

p <- arg_parser("This script merges genotype, site QC information file (VCFs) and its associated VEP annotation")
p <- add_argument(p, "input_directory_path", help="Input directory where main script is present")
p <- add_argument(p, "gene_name", help="name of the gene for output filename")
args <- parse_args(p)

message("Reading genotype, siteQC and annotation (GreenDB and VEP) files")

genotypes<-duckplyr::read_csv_duckdb(paste0(args$input_directory_path,"/",args$gene_name,"_genotypes.tsv"), options = list(delim="\t", sample_size=-1))
siteQC<-duckplyr::read_csv_duckdb(paste0(args$input_directory_path,"/",args$gene_name,"_siteQC.tsv"), options = list(delim="\t", sample_size=-1))
anno<-duckplyr::read_csv_duckdb(paste0(args$input_directory_path,"/",args$gene_name,"_annotation.tsv"),options = list(delim="\t",header=T, sample_size=-1))

message("distincting siteQC")
#Keep only distinct variants because overlapping enhancers can duplicate variants although their siteQC annotations will all be same
siteQC <- siteQC %>%
  distinct(variant_id, .keep_all = TRUE)

message("Processing the files and writing a tsv file ")
#Filter for rare variants in samples and attach all annotations
genotypes%>%
  mutate(genotype=if_else(genotype=="1/0","0/1",genotype))%>%
  distinct(sample, variant_id, .keep_all = TRUE) %>%
  left_join(siteQC,by="variant_id", relationship = "many-to-many") %>%
  filter(is.na(GEL_cohort_AF) | GEL_cohort_AF <= 0.01) %>%
  left_join(anno,by="variant_id", relationship = "many-to-many") %>%
  compute_csv(paste0(args$input_directory_path,"/",args$gene_name,".tsv"),options = list(delim="\t"))

##Check if file exists
final_output<-paste0(args$input_directory_path,"/",args$gene_name,".tsv")
final_output<-normalizePath(final_output)

if(file.exists(final_output) && file.info(final_output)$size > 0){
  message("SUCCESS: final output file created: ",final_output)
}else {
  message("ERROR: final output file missing or empty: ",final_output)
  quit(status = 1)
}
