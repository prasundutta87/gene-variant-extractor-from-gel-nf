#!/usr/bin/env Rscript
# join_additional_anno_with_genes.R
# Annotates the merged rare-variant table with participant metadata, phenotype
# data, variant classification flags, and per-variant carrier counts.
#
# Args (passed by Nextflow process, sourced from params):
#   1. gene_name              : HGNC gene symbol
#   2. sample_list            : AggV3 sample list CSV (duplicate platekey info)
#   3. sequencing_report_100k : GEL 100k sequencing report TSV
#   4. sequencing_report_gms  : GEL GMS sequencing report TSV
#   5. rd_phenotype           : RD v19 phenotype data TSV
#   6. gms_phenotype          : GMS v4 phenotype data TSV
#
# Input (in working directory, produced by combine_variant_using_duckplyr_with_siteQC.R):
#   ${gene}.tsv
#
# Output:
#   ${gene}/${gene}_for_review.tsv : annotated rare-variant table with carrier counts
#
# Notes:
#   - Duplicate platekeys: only the most recently delivered platekey per
#     participant is retained; all non-duplicated platekeys are kept as-is.
#   - affection_combined reconciles RD and GMS affection status; conflicting
#     sources are excluded from carrier counting.
#   - IMPACT filtering is intentionally absent — applied in gene-specific scripts.
#   - Exits with status 1 if output is missing or empty so Nextflow can detect failure.

library(tidyverse)
library(data.table)
library(gt)
library(ggpubr)
options(scipen=999)

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
  stop("Usage: join_additional_anno_with_genes.R <gene_name> <sample_list> ",
       "<sequencing_report_100k> <sequencing_report_gms> ",
       "<rd_phenotype> <gms_phenotype>")
}

gene_name             <- args[1]
sample_list_path      <- args[2]
seq_report_100k_path  <- args[3]
seq_report_gms_path   <- args[4]
input_RD_phenotype    <- args[5]
input_GMS_phenotype   <- args[6]

###########################################################################

##GEL AggV3 tells that-
#duplicate_of: We have 18 duplicates (n=36). This column provides the duplicated platekey of the duplicate pairing.
#If row is NA, no known duplicate has been identified at this stage.

participant_metadata<-fread(sample_list_path)

##these are 36 samples
duplicated_sample_list<-participant_metadata%>%filter(!is.na(duplicate_of))

##get sequencing reports from 100k and GMS data to identify the most recent sample out of a duplicate sample set
sequencing_report_100k<-fread(seq_report_100k_path)
sequencing_report_GMS<-fread(seq_report_gms_path)

sequencing_report_combined <- bind_rows(
  sequencing_report_100k %>%
    rename(platekey = plate_key) %>%
    select(platekey, delivery_date) %>%
    mutate(delivery_date = as.Date(delivery_date)),
  sequencing_report_GMS %>%
    select(platekey, delivery_date) %>%
    mutate(delivery_date = as.Date(delivery_date))
) %>%
  distinct(platekey, .keep_all = TRUE)

# For duplicate participants, identify the preferred (most recently delivered) platekey
preferred_platekeys <- duplicated_sample_list %>%
  left_join(sequencing_report_combined, by = "platekey") %>%
  group_by(pmin(platekey, duplicate_of)) %>%
  arrange(desc(delivery_date)) %>%
  slice(1) %>%
  ungroup() %>%
  pull(platekey)
##n=30 cases where unique platekeys match a single unique participant ID (GMS)
#unique_platekeys_matching_unique_participant_ID<-participant_metadata%>%count(participant_id)%>%
  #filter(n>1)

###########################################################################

# Input phenotype files
RD_phenotype_data<-fread(input_RD_phenotype, na.strings = c("", "NA", "."))
RD_phenotype_data<-RD_phenotype_data%>%mutate(participant_id=as.character(participant_id))
GMS_phenotype_data<-fread(input_GMS_phenotype, na.strings = c("", "NA", "."))
GMS_phenotype_data<-GMS_phenotype_data%>%mutate(participant_id=as.character(participant_id))

annotate_with_participant_metadata <- function(gene_name) {
  result <- fread(paste0(gene_name, ".tsv"),na.strings = c("", "NA", "."))%>%
    rename(platekey = sample) %>%
    mutate(
      DS = pmax(
        SpliceAI_pred_DS_AG,
        SpliceAI_pred_DS_DG,
        SpliceAI_pred_DS_AL,
        SpliceAI_pred_DS_DL
      )
    ) %>%
  mutate(
    DS_range = case_when(
      DS >= 0.2 & DS <= 0.5 ~ "0.2 to 0.5",
      DS > 0.5 & DS <= 0.8  ~ "0.5 to 0.8",
      DS > 0.8 & DS <= 1    ~ "0.8 to 1",
      TRUE ~ NA
    ),
    DS_source = case_when(
      is.na(DS) | DS == 0       ~ NA_character_,
      DS == SpliceAI_pred_DS_AG ~ "Acceptor Gain",
      DS == SpliceAI_pred_DS_DG ~ "Donor Gain",
      DS == SpliceAI_pred_DS_AL ~ "Acceptor Loss",
      DS == SpliceAI_pred_DS_DL ~ "Donor Loss",
      TRUE ~ NA_character_
    ),
    zygosity = case_when(
      genotype %in% c("1/1", "1") ~ "HOM/HEMI",
      genotype == "0/1"           ~ "HET"
    )
  ) %>%
  mutate(
    HGVSc = na_if(HGVSc, "-"),
    HGVSc_empty = is.na(HGVSc),
    SNV_with_valid_HGVSc = !HGVSc_empty &
      grepl("\\+|\\-", HGVSc) & !grepl("ins|del|dup", HGVSc),
    distance_from_nearest_exon = as.numeric(stringr::str_match(HGVSc, "[\\+\\-](\\d+)")[, 2]),
    deep_intronic = !is.na(distance_from_nearest_exon) &
      distance_from_nearest_exon > 20,
    HGVSc_annotation_category = case_when(
      is.na(HGVSc)                      ~ "Missing HGVSc",
      grepl("ins|del|dup", HGVSc)       ~ "Indel/dup",
      !grepl("\\+|\\-", HGVSc)          ~ "Not intronic",
      is.na(distance_from_nearest_exon) ~ "Cannot extract distance",
      .default                          = "Valid intronic SNV"
    )
  ) %>%
  left_join(participant_metadata,by = "platekey")%>%
  relocate(participant_id,type,study_source,programme,duplicate_of) %>%
  select(
    -starts_with("gnomAD"),
    gnomADg,
    gnomADg_AF_exomes,
    gnomADg_AF_exomes_XX,
    gnomADg_AF_exomes_XY,
    gnomADg_AF_genomes,
    gnomADg_AF_genomes_XX,
    gnomADg_AF_genomes_XY,
    gnomADg_AF_joint,
    gnomADg_AF_joint_XX,
    gnomADg_AF_joint_XY
  ) %>%
  left_join(RD_phenotype_data,by="participant_id")%>%
  left_join(GMS_phenotype_data,by="participant_id")
  return(result)
}

gene_data<-annotate_with_participant_metadata(gene_name)

#' normalise_affection_status
#'
#' Integrates two sources of affection status (RD and GMS) into a single
#' combined label per participant. Agreement between sources yields
#' "affected" or "unaffected". Disagreement yields "conflicting" which
#' is excluded from downstream segregation scoring. Single-source
#' participants retain that source's label.
#'
#' @param data Input dataframe containing RD_affection_status and
#'             GMS_disease_status columns
#'
#' @return Original dataframe with added affection_combined column
#'         (rd_status and gms_status are dropped as intermediate columns)
#'  - "affected"    : both sources agree affected, or single source affected
#'  - "unaffected"  : both sources agree unaffected, or single source unaffected
#'  - "conflicting" : sources disagree — excluded from scoring
#'  - NA            : no affection data available from either source

normalise_affection_status <- function(df) {
  df %>%
    mutate(
      rd_status  = tolower(RD_affection_status),
      gms_status = tolower(GMS_disease_status),
      affection_combined = case_when(
        rd_status == "affected"   & gms_status == "affected"   ~ "affected",
        rd_status == "unaffected" & gms_status == "unaffected" ~ "unaffected",
        rd_status == "affected"   & is.na(gms_status)          ~ "affected",
        rd_status == "unaffected" & is.na(gms_status)          ~ "unaffected",
        is.na(rd_status)          & gms_status == "affected"   ~ "affected",
        is.na(rd_status)          & gms_status == "unaffected" ~ "unaffected",
        rd_status == "affected"   & gms_status == "unaffected" ~ "conflicting",
        rd_status == "unaffected" & gms_status == "affected"   ~ "conflicting",
        TRUE                                                    ~ NA_character_
      )
    ) %>%
    select(-rd_status, -gms_status)
}

#' classify_variants_for_review
#'
#' Classifies variants by AF, ClinVar status, and call quality.
#' Designed to be universal across genes — no IMPACT filtering applied here.
#' Gene-specific filtering (IMPACT, SpliceAI, inheritance) should be applied
#' downstream in gene-specific scripts.
#'
#' @param data Input dataframe — one row per participant per variant
#' @param very_rare_af AF threshold for "very rare" (default: 0.00001)
#' @param rare_af AF threshold for "rare" (default: 0.0001)
#'
#' @return Original dataframe with added columns:
#'  - consequence_clean, clinvar_clean    : normalised text fields
#'  - gnomad_af_value, gel_af_value       : numeric AF values
#'  - gnomad_scenario, gel_scenario       : descriptive AF categories
#'  - gnomad_compatible, gel_compatible   : AF compatibility flags
#'  - is_clinvar_benign                   : ClinVar benign flag
#'  - keep_for_review                     : FALSE if benign or too common
#'  - review_note                         : readable explanation
#'  - passes_quality, quality_note        : call quality flag and note
#'  - clinvar_pathogenic_failed_quality   : pathogenic variants failing QC

classify_variants_for_review <- function(data, very_rare_af = 0.00001, rare_af = 0.0001) {
  data %>%
    mutate(
      consequence_clean = str_to_lower(str_replace_all(Consequence, "\\s+", "")),
      clinvar_clean     = str_to_lower(str_replace_all(ClinVar_CLNSIG, "_", " ")),
      gnomad_af_value   = suppressWarnings(as.numeric(gnomADg_AF_joint)),
      gel_af_value      = suppressWarnings(as.numeric(GEL_cohort_AF)),
      is_clinvar_benign = !is.na(clinvar_clean) & str_detect(clinvar_clean, "benign"),
      gnomad_scenario = case_when(
        is.na(gnomad_af_value) | gnomad_af_value == 0 ~ "absent_or_not_in_gnomad",
        gnomad_af_value <= very_rare_af                ~ "very_rare",
        gnomad_af_value <= rare_af                     ~ "rare",
        gnomad_af_value > rare_af                      ~ "too_common",
        TRUE                                           ~ "unknown"
      ),
      gnomad_compatible = gnomad_scenario %in% c("absent_or_not_in_gnomad", "very_rare", "rare"),
      gel_scenario = case_when(
        is.na(gel_af_value) | gel_af_value == 0 ~ "absent_in_gel",
        gel_af_value <= very_rare_af             ~ "very_rare_in_gel",
        gel_af_value <= rare_af                  ~ "rare_in_gel",
        gel_af_value > rare_af                   ~ "common_in_gel",
        TRUE                                     ~ "unknown"
      ),
      gel_compatible  = gel_scenario %in% c("absent_in_gel", "very_rare_in_gel", "rare_in_gel"),
      keep_for_review = case_when(
        is_clinvar_benign  ~ FALSE,
        !gnomad_compatible ~ FALSE,
        !gel_compatible    ~ FALSE,
        TRUE               ~ TRUE
      ),
      review_note = case_when(
        is_clinvar_benign  ~ paste0("excluded: ClinVar=", ClinVar_CLNSIG),
        !gnomad_compatible ~ paste0("excluded: gnomAD AF=", round(gnomad_af_value, 6), " | ", gnomad_scenario),
        !gel_compatible    ~ paste0("excluded: GEL cohort AF=", round(gel_af_value, 6), " | ", gel_scenario),
        keep_for_review    ~ paste0("kept: ", IMPACT, " | ", Consequence, " | gnomAD=", gnomad_scenario, " | GEL=", gel_scenario),
        TRUE               ~ paste0("not kept: ", Consequence)
      ),
      passes_quality = case_when(
        FILTER           != "PASS" ~ FALSE,
        medianDP         < 10      ~ FALSE,
        medianGQ         < 20      ~ FALSE,
        missingness_rate > 0.05   ~ FALSE,
        TRUE                       ~ TRUE
      ),
      quality_note = case_when(
        FILTER           != "PASS" ~ paste0("failed: variant caller | FILTER = ", FILTER),
        medianDP         < 10      ~ paste0("failed: low depth | DP = ", medianDP),
        medianGQ         < 20      ~ paste0("failed: low GQ | GQ = ", medianGQ),
        missingness_rate > 0.05   ~ paste0("failed: high missingness | missingness_rate = ", round(missingness_rate, 4)),
        TRUE                       ~ paste0("passed: DP = ", medianDP, " | GQ = ", medianGQ, " | missingness = ", round(missingness_rate, 4))
      )
    ) %>%
    mutate(
      clinvar_pathogenic_failed_quality = str_detect(clinvar_clean, "pathogenic") &
        !str_detect(clinvar_clean, "conflicting") &
        !is_clinvar_benign &
        !passes_quality
    )
}

normalise_metadata <- function(data) {
  data %>%
    mutate(
      # -- Sex ----------------------------------------
      sex = case_when(
        study_source == "GMS"    ~ str_to_lower(GMS_administrative_gender),
        study_source == "MP100K" ~ str_to_lower(RD_participant_phenotypic_sex),
        TRUE ~ NA_character_
      ),
      # -- Participant type ----------------------------
      participant_type = case_when(
        study_source == "MP100K"                  ~ RD_participant_type,
        GMS_participant_type == "FullSibling"     ~ "Full Sibling",
        GMS_participant_type == "TwinsMonozygous" ~ "Twins Monozygous",
        study_source == "GMS"                     ~ GMS_participant_type,
        TRUE ~ NA_character_
      ),
      # -- Case solved --------------------------------
      case_solved = case_when(
        study_source == "GMS"    ~ GMS_exit_questionnaire_case_solved_family,
        study_source == "MP100K" ~ RD_gmc_exit_case_solved_family,
        TRUE ~ NA_character_
      ),
      # -- HPO terms ----------------------------------
      hpo_terms = case_when(
        study_source == "GMS"    ~ GMS_normalised_hpo_term,
        study_source == "MP100K" ~ RD_normalised_hpo_term,
        TRUE ~ NA_character_
      ),
      # -- Year of birth / Age ------------------------
      year_of_birth = case_when(
        study_source == "MP100K" ~ as.numeric(RD_yob),
        study_source == "GMS"    ~ as.numeric(GMS_participant_year_of_birth),
        TRUE ~ NA_real_
      ),
      age = as.integer(format(Sys.Date(), "%Y")) - year_of_birth,
      age_reference_year = as.integer(format(Sys.Date(), "%Y"))
    ) %>%
    select(
      -GMS_administrative_gender,
      -RD_participant_phenotypic_sex,
      -RD_participant_type,
      -GMS_participant_type,
      -GMS_exit_questionnaire_case_solved_family,
      -RD_gmc_exit_case_solved_family,
      -GMS_normalised_hpo_term,
      -RD_normalised_hpo_term,
      -RD_yob,
      -GMS_participant_year_of_birth
    )
}

# FUNCTION: count_carriers_by_affection
#
# Counts carriers per variant broken down by zygosity and affection status.
# Filters to MANE_SELECT, quality-passed, rare_disease programme rows only.
# Run after normalise_affection_status() and classify_variants_for_review().
#
# @param df    Dataframe with affection_combined, zygosity, passes_quality columns
# @return      Summarised dataframe with one row per variant containing:
#              n_affected_het/hom, affected_het/hom_ids, n_unaffected_het/hom,
#              n_conflicting, n_missing, n_total

count_carriers_by_affection <- function(df) {
  df %>%
    filter(!is.na(MANE_SELECT)) %>%
    filter(passes_quality == TRUE) %>%
    filter(programme == "rare_disease") %>%
    #Safety net: ensure one row per participant per variant
    distinct(variant_id, participant_id, .keep_all = TRUE) %>%
    group_by(variant_id) %>%
    summarise(
      n_affected_het   = sum(affection_combined == "affected"   & zygosity == "HET",      na.rm = TRUE),
      affected_het_ids = paste(participant_id[affection_combined == "affected"   & zygosity == "HET"],     collapse = ";"),
      n_affected_hom   = sum(affection_combined == "affected"   & zygosity == "HOM/HEMI", na.rm = TRUE),
      affected_hom_ids = paste(participant_id[affection_combined == "affected"   & zygosity == "HOM/HEMI"], collapse = ";"),
      n_unaffected_het = sum(affection_combined == "unaffected" & zygosity == "HET",      na.rm = TRUE),
      n_unaffected_hom = sum(affection_combined == "unaffected" & zygosity == "HOM/HEMI", na.rm = TRUE),
      n_conflicting    = sum(affection_combined == "conflicting", na.rm = TRUE),
      n_missing        = sum(is.na(affection_combined)),
      n_total          = n(),
      .groups = "drop"
    )
}

#Add additional phenotypic annotation and calculate various flags to be used downstream for filtering
gene_annotated <- gene_data %>%
  filter(
    !platekey %in% duplicated_sample_list$platekey |  # not in duplicate list - keep all
      platekey %in% preferred_platekeys,              # in duplicate list - keep only most recently delivered
    programme == "rare_disease"                       ##restrict to rare diseases participants
  ) %>%
  normalise_affection_status() %>%
  normalise_metadata() %>%
  classify_variants_for_review()

# Step 3 - count carriers on quality-passed calls keeping rare diseases participants in mind only
gene_affection_counts <- gene_annotated %>%
  count_carriers_by_affection()

# Step 4 - join counts back
gene_for_review <- gene_annotated %>%
  left_join(gene_affection_counts, by = "variant_id")

dir.create(gene_name, showWarnings = FALSE)
fwrite(gene_for_review,
       paste0(gene_name, "/", tolower(gene_name), "_for_review.tsv.gz"),
       sep = "\t", row.names = FALSE, quote = FALSE, compress = "gzip")

# ── Output validation ─────────────────────────────────────────────────────────
output_path <- paste0(gene_name, "/", tolower(gene_name), "_for_review.tsv.gz")
if (file.exists(output_path) && file.info(output_path)$size > 0) {
  message("SUCCESS: output written to ", output_path)
} else {
  message("ERROR: output file missing or empty: ", output_path)
  quit(status = 1)
}

############################################SCRIPT END############################################
