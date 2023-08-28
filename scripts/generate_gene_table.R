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
        input=list(sce="data/sce/merged.genes.full_annot.Rds"),
        output=list(cpc="/fscratch/fanslerm/merged_gene_cpc_pointestimates.tsv.gz",
                    tpm="/fscratch/fanslerm/merged_gene_tpm_pointestimates.tsv.gz"),
        params=list())
}

################################################################################
## Load Data
################################################################################

## Load UTR Counts
sce <- readRDS(snakemake@input$sce)

################################################################################
## DESIGN MATRICES 
################################################################################

## M: (cells) x (cell types)
M_cells_celltypes <- colData(sce)[,c('tissue', 'cell_type', 'age')] %>%
    as_tibble %>%
    mutate(tissue_celltype_age=paste(
               str_replace_all(tissue, "_", " "),
               cell_type,
               age, sep=', ')) %>%
    pull(tissue_celltype_age) %>%
    fac2sparse %>%
    t

################################################################################
## Total Gene Expression 
################################################################################

cts_genes_celltypes <- assay(sce, 'normcounts') %*% M_cells_celltypes

## n cells per celltype
ncells_celltypes <- colSums(M_cells_celltypes)

cpc_genes_celltypes <- cts_genes_celltypes %>%
    { . %*% Diagonal(ncol(.), 1/ncells_celltypes) } %>%
    `colnames<-`(colnames(cts_genes_celltypes))

tpm_genes_celltypes <- cts_genes_celltypes %>%
    { . %*% Diagonal(ncol(.), 1e6/colSums(.)) } %>%
    `colnames<-`(colnames(cts_genes_celltypes))

df_genes <- rowData(sce) %>%
    as_tibble

cpc_genes_celltypes %>%
    as.matrix %>% as.data.frame %>%
    rownames_to_column('gene_id') %>%
    right_join(x=df_genes, by='gene_id') %>%
    arrange(gene_name) %>%
    write_tsv(snakemake@output$cpc) ## Export

tpm_genes_celltypes %>%
    as.matrix %>% as.data.frame %>%
    rownames_to_column('gene_id') %>%
    right_join(x=df_genes, by='gene_id') %>%
    arrange(gene_name) %>%
    write_tsv(snakemake@output$tpm) ## Export
