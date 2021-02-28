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
## Mock Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(sce="data/sce/merged.genes.raw.Rds",
                   utrs="data/utrs/utrome_genes_annotation.tsv",
                   size_factors="data/scran/merged.size_factors.tsv.gz"),
        output=list(sce="/fscratch/fanslerm/merged.genes.full_annot.Rds"),
        params=list())
}

################################################################################
## Load Data
################################################################################

sce <- readRDS(snakemake@input$sce)

df_utrs <- read_tsv(snakemake@input$utrs, col_types='ccicilllcc') %>%
    set_rownames(.$gene_id) %>%
    DataFrame()

df_sizeFactors <- read_tsv(snakemake@input$size_factors, col_types='cdd', skip=1,
                           col_names=c('cell_id', 'atlas.size_factor', 'atlas.capture_efficiency'))

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
