#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(SingleCellExperiment)
library(S4Vectors)
library(tidyverse)


################################################################################
## Load Arguments
################################################################################

if (interactive()) {
    args <- c("data/sce/TabulaMuris.genes.annot.rds",
              "data/utrs/genes-utr-metadata-lengths.tsv",
              "data/scran/merged.size_factors.tsv.gz",
              "/scratch/fanslerm/TM.genes.full_annot.rds")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 4) {
        stop("Incorrect number of arguments!\nUsage:\n> annotate_genes_sce.R <sceRDS> <utrsFile> <sizesFile> <outFile>\n")
    }
}
arg.sceRDS <- args[1]
arg.utrFile <- args[2]
arg.sizeFactors <- args[3]
arg.outFile <- args[4]


################################################################################
## Load Data
################################################################################

sce <- readRDS(arg.sceRDS)

df.utrs <- read_tsv(arg.utrFile, col_types='_ccidcdi')

df.sizeFactors <- read_tsv(arg.sizeFactors, col_types='cdd') %>%
    rename_at(vars(-c('cell')), ~ str_c(., '.merged'))


################################################################################
## Reannotate Genes
################################################################################

## Save previous labels for validation
old.rownames <- rownames(sce)
old.colnames <- colnames(sce)

rowData(sce) %<>%
    as.data.frame() %>%
    rownames_to_column() %>%
    select(-one_of(c('utr.count.mca', 'utr.count.tm', 'utr.type.mca', 'utr.type.tm', 'utr.type.matches'))) %>%
    left_join(df.utrs, by='gene_id') %>%
    column_to_rownames() %>%
    DataFrame()

colData(sce) %<>%
    as.data.frame() %>%
    rownames_to_column() %>%
    left_join(df.sizeFactors, by="cell") %>%
    column_to_rownames() %>%
    DataFrame()

## Ensure we are not scrambling genes or cells
stopifnot(all(rownames(sce) == old.rownames))
stopifnot(all(colnames(sce) == old.colnames))


################################################################################
## Export
################################################################################

saveRDS(sce, arg.outFile)
