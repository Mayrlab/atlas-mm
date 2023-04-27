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
        input=list(utrs="data/utrs/utrome_txs_annotation.tsv"),
        output=list(tsv="/fscratch/fanslerm/utrome_genes_annotation.tsv",
                    rds="/fscratch/fanslerm/utrome_genes_annotation.Rds"),
        params=list())
}

################################################################################
## Load Data
################################################################################

df_utrs <- read_tsv(snakemake@input$utrs, col_types='cccccilllildidiicill')

################################################################################
## Collapse to Genes
################################################################################

utr_lengths_str <- function (utr_length, idx) {
    if (any(idx)) {
        str_c(utr_length[idx], collapse=";")
    } else { "" }
}


df_genes_final <- df_utrs %>%
    group_by(gene_id, gene_name, atlas.ncelltypes_gene, atlas.utr_type,
             atlas.n_utrs_no_ipa, is_blacklisted) %>%
    arrange(utr_length) %>%
    summarize(is_consistent=all(is_consistent), has_ipa=any(is_ipa),
              utr_lengths_tandem=utr_lengths_str(utr_length, !is_ipa),
              utr_lengths_ipa=utr_lengths_str(utr_length, is_ipa),
              .groups='drop')

################################################################################
## Export
################################################################################

write_tsv(df_genes_final, snakemake@output$tsv)

df_genes_final %>%
    DataFrame(row.names=.$gene_id) %>%
    saveRDS(snakemake@output$rds)
