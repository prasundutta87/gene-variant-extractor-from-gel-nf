#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(argparser)))
suppressMessages(suppressWarnings(library(duckplyr)))

p <- arg_parser("This script merges genotype, site QC information file (VCFs) and its associated VEP annotation")

#Add Positional arguments
p <- add_argument(p, "input_directory_path", help="Input directory where main script is present")
p <- add_argument(p, "gene_name", help="name of the gene for output filename")
args <- parse_args(p)

genotypes<-duckplyr::read_csv_duckdb(paste0(args$input_directory_path,"/",args$gene_name,"_genotypes.tsv"), options = list(delim="\t", sample_size=-1))
siteQC<-duckplyr::read_csv_duckdb(paste0(args$input_directory_path,"/",args$gene_name,"_siteQC.tsv"), options = list(delim="\t", sample_size=-1))
anno<-duckplyr::read_csv_duckdb(paste0(args$input_directory_path,"/",args$gene_name,"_annotation.tsv"),options = list(delim="\t",header=T, sample_size=-1))

genotypes%>%
  mutate(genotype=if_else(genotype=="1/0","0/1",genotype))%>%
  left_join(siteQC,by="variant_id", relationship = "many-to-many") %>%
  left_join(anno,by="variant_id", relationship = "many-to-many") %>% 
  compute_csv(paste0(args$input_directory_path,"/",args$gene_name,".tsv"),options = list(delim="\t"))

##Check if file exists

final_output<-paste0(args$input_directory_path,"/",args$gene_name,".tsv")

final_output<-normalizePath(final_output)

if(file.exists(final_output) && file.info(final_output)$size > 5000){
 message("SUCCESS: final output file created: ",final_output)
}else {
  message("ERROR: final output file missing or empty: ",final_output)
  }