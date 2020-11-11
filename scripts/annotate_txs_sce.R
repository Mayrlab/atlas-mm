#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(SingleCellExperiment)
library(S4Vectors)
library(tidyverse)
library(magrittr)

################################################################################
## Mock Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(sce="/data/mayrc/data/tmuris/sce/tmuris.utrome.txs.Rds",
                   utrs="data/utrs/txs_utr_metadata_lengths.tsv",
                   size_factors="data/scran/merged.size_factors.tsv.gz"),
        output=list(sce="/fscratch/fanslerm/tmuris.utrome.txs.full_annot.Rds"),
        params=list())
}

################################################################################
## Load Data
################################################################################

sce <- readRDS(snakemake@input$sce)

df.utrs <- read_tsv(snakemake@input$utrs, col_types='c____di__l__d___ii___i_ccil') %>%
    dplyr::rename(utr.ncelltypes.no_ipa=utr.ncelltypes.pct10.no_ipa,
                  utr.count.no_ipa=utr.count.pct10.no_ipa,
                  utr.type.no_ipa=utr.type.pct10.no_ipa,
                  improper_utr_length=improper)

df.sizeFactors <- read_tsv(snakemake@input$size_factors, col_types='cdd') %>%
    rename_at(vars(-c('cell_id')), ~ str_c(., '.merged'))


################################################################################
## Reannotate Transcripts
################################################################################

## Save previous labels for validation
old.rownames <- rownames(sce)
old.colnames <- colnames(sce)

rowData(sce) %<>%
    as.data.frame() %>%
    rownames_to_column() %>%
    left_join(df.utrs, by='transcript_id') %>%
    column_to_rownames() %>%
    DataFrame()

colData(sce) %<>%
    as.data.frame() %>%
    mutate(cell_id=rownames(.)) %>%
    left_join(df.sizeFactors, by="cell_id") %>%
    set_rownames(.$cell_id) %>%
    DataFrame()

## Ensure we are not scrambling genes or cells
stopifnot(all(rownames(sce) == old.rownames))
stopifnot(all(colnames(sce) == old.colnames))

################################################################################
## Export
################################################################################

saveRDS(sce, snakemake@output$sce)
