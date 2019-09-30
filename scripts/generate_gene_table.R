#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(SingleCellExperiment)
library(tidyverse)
library(Matrix)


################################################################################
## Load Arguments
################################################################################

if (interactive()) {
    args <- c("data/sce/merged.genes.full_annot.rds",
              "data/blacklist.utrome.txt",
              "/scratch/fanslerm/TM.gene.tsv")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 3) {
        stop("Incorrect number of arguments!\nUsage:\n> generate_gene_table.R <sceRDS> <blacklist> <outFile>\n")
    }
}
arg.sceRDS    <- args[1]
arg.blacklist <- args[2]
arg.outFile   <- args[3]

################################################################################
## Load Data
################################################################################

## Load UTR Counts
sce <- readRDS(arg.sceRDS)
sce <- sce[, !is.na(sce$size_factor.merged)]

## Known blacklist
genes.blacklist <- read_lines(arg.blacklist)

################################################################################
## DESIGN MATRICES 
################################################################################

## M: (cells) x (cell types)
M.celltypes <- colData(sce)[,c('tissue', 'cell_type', 'age')] %>%
    as.data.frame() %>%
    mutate(tissue_celltype_age=paste(
               str_replace_all(tissue, "_", " "),
               cell_type,
               age, sep=', ')) %>%
    pull(tissue_celltype_age) %>%
    fac2sparse() %>%
    t()

################################################################################
## Total Gene Expression 
################################################################################

idx.genes <- rowData(sce) %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(expressed.gene, !(gene_symbol %in% genes.blacklist)) %>%
    pull(rowname)

cts.clusters.all <- counts(sce[idx.genes,]) %*% (M.celltypes / sce$size_factor.merged)

mu.gene.clusters.all <- t(t(cts.clusters.all) / colSums(M.celltypes))

df.genes <- rowData(sce) %>%
    as.data.frame() %>%
    select('gene_symbol', 'gene_id', 'utr_count')

mu.gene.clusters.all %>%
    as.matrix() %>% as.data.frame() %>%
    rownames_to_column('gene_id') %>%
    left_join(df.genes, by='gene_id') %>%
    select(gene_symbol, gene_id, utr_count, everything()) %>%
    arrange(gene_symbol) %>%
    write_tsv(arg.outFile) ## Export
