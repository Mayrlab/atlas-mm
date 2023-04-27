#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(tidyverse)
library(magrittr)
library(S4Vectors)

################################################################################
## Mock Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(utrs="data/utrs/txs_utr_metadata_lengths.tsv"),
        output=list(tsv="/fscratch/fanslerm/utrome_txs_annotation.tsv",
                    rds="/fscratch/fanslerm/utrome_txs_annotation.Rds"),
        params=list())
}

################################################################################
## Load Data
################################################################################

df_utrs <- read_tsv(snakemake@input$utrs, col_types='ccc__di__l__d___ii___i_cccill') %>%
    dplyr::rename(atlas.pct_utr_total=utr.pct.total,
                  atlas.rank_utr_total=utr.rank.total,
                  atlas.pct_utr_no_ipa=utr.pct.no_ipa,
                  atlas.ncelltypes_utr=utr.ncelltypes.pct10.no_ipa,
                  atlas.ncelltypes_gene=gene.ncelltypes.cells50.no_ipa,
                  atlas.n_utrs_no_ipa=utr.count.pct10.no_ipa,
                  atlas.utr_type=utr.type.pct10.no_ipa,
                  is_improper_utr_length=improper) %>%
    mutate(utr_position=as.integer(str_extract(utr_name, "[0-9]+$"))) 

################################################################################
## Annotation Distal and LU
################################################################################

df_lus <- df_utrs %>%
    filter(!is_ipa, (atlas.ncelltypes_utr > 0) | (atlas.pct_utr_no_ipa >= 0.1)) %>%
    group_by(gene_id) %>%
    mutate(is_distal=utr_position == max(utr_position),
           ## if length not computed, fall back to distal
           is_lu=ifelse(is.na(utr_length), is_distal, utr_length == max(utr_length, na.rm=TRUE))) %>%
    ungroup()

txs_lu <- filter(df_lus, is_lu) %$% transcript_id
txs_distal <- filter(df_lus, is_distal) %$% transcript_id

df_utrs_final <- df_utrs %>%
    mutate(is_lu=transcript_id %in% txs_lu,
           is_distal=transcript_id %in% txs_distal) %>%
    group_by(gene_id) %>%
    mutate(is_distal=if (any(is_distal)) is_distal else utr_position == max(utr_position),
           is_lu=if (any(is_lu)) is_lu else is_distal,
           is_consistent=is_distal == is_lu) %>%
    ungroup() %>%
    select("transcript_id", "transcript_gencode", "gene_id", "gene_name", "utr_name",
           "utr_position", "is_ipa", "is_lu", "is_distal",
           "utr_length", "is_improper_utr_length", 
           "atlas.pct_utr_total", "atlas.rank_utr_total",
           "atlas.pct_utr_no_ipa", "atlas.ncelltypes_utr", 
           "atlas.ncelltypes_gene", "atlas.utr_type", "atlas.n_utrs_no_ipa", 
           "is_blacklisted", "is_consistent")


################################################################################
## Export
################################################################################

write_tsv(df_utrs_final, snakemake@output$tsv)

df_utrs_final %>%
    DataFrame(row.names=.$transcript_id) %>%
    saveRDS(snakemake@output$rds)
