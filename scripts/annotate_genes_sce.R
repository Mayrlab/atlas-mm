#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(SingleCellExperiment)
library(scater)
library(S4Vectors)
library(tidyverse)
library(magrittr)

################################################################################
## Mock data
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(
            ##sce="/data/mayrc/data/tmuris/sce/tmuris.utrome.genes.Rds",
            ##sce="/data/mayrc/data/mm-hspc-44k/sce/hspcs.utrome.genes.Rds",
            sce="/data/mayrc/fanslerm/scutr-quant/data/sce/ximerakis19.utrome.txs.Rds",
            utrs="data/utrs/utrome_genes_annotation.Rds",
            size_factors="data/scran/merged.size_factors.tsv.gz"),
        output=list(sce="/fscratch/fanslerm/tmuris.utrome.txs.full_annot.Rds"),
        params=list())
}

################################################################################
## Load Data
################################################################################

sce <- readRDS(snakemake@input$sce)

df_utrs <- readRDS(snakemake@input$utrs)

df_sizeFactors <- read_tsv(snakemake@input$size_factors, col_types='cdd', skip=1,
                           col_names=c('cell_id', 'atlas.size_factor', 'atlas.capture_efficiency'))

## Experiement-Specific Adjustments
## brain SCE includes unlabeled cells, let's remove
if ('cell_type' %in% colnames(colData(sce))) {
    sce %<>% `[`(, !is.na(.$cell_type))
}

## HSPCs have cell_id's misordered
if (!any(colnames(sce) %in% df_sizeFactors$cell_id)) {
    colData(sce) %<>%
        as_tibble %>%
        mutate(cell_id=str_replace(cell_id, '^([ACGT]{16})_(.*)$', '\\2_\\1')) %>%
        set_rownames(.$cell_id) %>%
        DataFrame()
}

################################################################################
## Reannotate Genes
################################################################################

## Save previous labels for validation
old.rownames <- rownames(sce)
old.colnames <- colnames(sce)

rowData(sce) <- df_utrs[old.rownames,]

colData(sce) %<>%
    as.data.frame() %>%
    mutate(cell_id=rownames(.)) %>%
    left_join(df_sizeFactors, by="cell_id") %>%
    set_rownames(.$cell_id) %>%
    DataFrame()

## Ensure we are not scrambling genes or cells
stopifnot(all(rownames(sce) == old.rownames))
stopifnot(all(colnames(sce) == old.colnames))

################################################################################
## Size Factors
################################################################################

sce %<>% `[`(,!is.na(.$atlas.size_factor))

sizeFactors(sce) <- sce$atlas.size_factor

assay(sce, 'normcounts') <- normalizeCounts(sce, log=FALSE)

################################################################################
## Export
################################################################################

saveRDS(sce, snakemake@output$sce)
