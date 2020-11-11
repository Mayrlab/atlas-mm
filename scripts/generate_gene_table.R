#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(SingleCellExperiment)
library(tidyverse)
library(Matrix)

################################################################################
## Mock Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(sce="data/sce/merged.genes.full_annot.Rds",
                   blacklist="data/blacklist.utrome.txt"),
        output=list(tsv="/fscratch/fanslerm/merged_gene_expr_pointestimates.tsv.gz"),
        params=list())
}

################################################################################
## Load Data
################################################################################

## Load UTR Counts
sce <- readRDS(snakemake@input$sce)
sce <- sce[, !is.na(sce$size_factor.merged)]

## Known blacklist
genes.blacklist <- read_lines(snakemake@input$blacklist)

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
    filter(!is.na(utr_type)) %>%
    pull(rowname)

cts.clusters.all <- counts(sce[idx.genes,]) %*% (M.celltypes / sce$size_factor.merged)

mu.gene.clusters.all <- t(t(cts.clusters.all) / colSums(M.celltypes))

df.genes <- rowData(sce) %>%
    as_tibble(rownames="gene_id") %>%
    select('gene_id', 'utr_type', 'utr_count')

mu.gene.clusters.all %>%
    as.matrix() %>% as.data.frame() %>%
    rownames_to_column('gene_id') %>%
    left_join(df.genes, by='gene_id') %>%
    select(gene_id, utr_type, utr_count, everything()) %>%
    arrange(gene_id) %>%
    write_tsv(snakemake@output$tsv) ## Export
