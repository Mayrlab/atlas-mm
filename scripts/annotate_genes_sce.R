#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(SingleCellExperiment)
library(S4Vectors)
library(tidyverse)
library(magrittr)

################################################################################
## Mock data
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(sce="/data/mayrc/data/tmuris/sce/tmuris.utrome.genes.Rds",
                   utrs="data/utrs/genes_utr_metadata_lengths.tsv",
                   size_factors="data/scran/merged.size_factors.tsv.gz"),
        output=list(sce="/fscratch/fanslerm/tmuris.utrome.txs.full_annot.Rds"),
        params=list())
}

################################################################################
## Load Data
################################################################################

sce <- readRDS(snakemake@input$sce)

df.utrs <- read_tsv(snakemake@input$utrs, col_types='_ccidcdi')

df.sizeFactors <- read_tsv(snakemake@input$size_factors, col_types='cdd') %>%
    rename_at(vars(-c('cell_id')), ~ str_c(., '.merged'))


################################################################################
## Reannotate Genes
################################################################################

## Save previous labels for validation
old.rownames <- rownames(sce)
old.colnames <- colnames(sce)

rowData(sce) %<>%
    as_tibble(rownames='gene_id') %>%
    left_join(df.utrs, by='gene_id') %>%
    column_to_rownames('gene_id') %>%
    DataFrame

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
