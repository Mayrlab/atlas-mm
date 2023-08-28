#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(SingleCellExperiment)
library(tidyverse)
library(magrittr)
library(S4Vectors)
library(Matrix)
library(matrixStats)

################################################################################
## Mock Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(sce="data/sce/merged.txs.full_annot.Rds",
                   genes="data/utrs/utrome_genes_annotation.tsv"),
        output=list(lui="/fscratch/fanslerm/merged_lui_expr_pointestimates.tsv",
                    n_cells="/fscratch/fanslerm/merged_ncells_expr.tsv"),
        params=list(min_cells=50))
}

################################################################################
## Load Data
################################################################################

## Load UTR Counts
sce <- readRDS(snakemake@input$sce)

df_genes <- read_tsv(snakemake@input$genes, col_types='cci_ill_c_')

################################################################################
## Identify Multiutrs
################################################################################

df_multiutrs <- rowData(sce) %>%
    as_tibble %>%
    filter(atlas.utr_type == 'multi',
           is_ipa == FALSE, is_blacklisted == FALSE,
           atlas.ncelltypes_gene >= 2) %>%
    arrange(gene_name, utr_position)

cat(sprintf("Considering **%d genes**.\n",
            df_multiutrs %>% pull(gene_id) %>% unique() %>% length()))

idx_utrs <- df_multiutrs %>% pull(transcript_id)

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

## M: (genes) x (utrs)
M_genes_txs <- rowData(sce[idx_utrs,])$gene_id %>% fac2sparse()
colnames(M_genes_txs) <- idx_utrs

## M: (distal UTRs) x (utrs)
M_distal_txs <- M_genes_txs %*% Diagonal(ncol(M_genes_txs), df_multiutrs$is_distal)

################################################################################
## COMPUTE LUI AND GENE COUNTS 
################################################################################

## transcripts counts
cts_txs_celltypes <- assay(sce[idx_utrs,], 'normcounts') %*% M_cells_celltypes

## n celltypes
ncells_celltypes <- colSums(M_cells_celltypes)

## mean gene counts
cpc_genes_celltypes <- (M_genes_txs %*% cts_txs_celltypes %*% Diagonal(length(ncells_celltypes), 1/ncells_celltypes)) %>%
    `colnames<-`(colnames(cts_txs_celltypes))

## number of cells in each cell type expressing a given gene
ncells_genes_celltypes <- ((M_genes_txs %*% counts(sce[idx_utrs,])) > 0) %*% M_cells_celltypes

## Point Estimates for LUI for all (gene) x (cell type)
lui_genes_celltypes <- (M_distal_txs %*% cts_txs_celltypes) / (M_genes_txs %*% cts_txs_celltypes)

## Mask all gene-cell combinations
lui_genes_celltypes[which(ncells_genes_celltypes < snakemake@params$min_cells)] <- NA

################################################################################
## Wide table 
################################################################################

df_lui <- lui_genes_celltypes %>%
    as.matrix %>% as.data.frame %>%
    rownames_to_column('gene_id') %>%
    select(gene_id, everything()) 


df_cpc <- cpc_genes_celltypes %>%
    as.matrix %>% as.data.frame %>%
    rownames_to_column('gene_id') %>%
    select(gene_id, everything())

df_full <- df_lui %>%
    left_join(df_cpc, by='gene_id', suffix=c(" LUI", " Mean Expression")) %>%
    right_join(x=df_genes, by='gene_id')

################################################################################
## Export
################################################################################

write_tsv(df_full, snakemake@output$lui)

ncells_genes_celltypes %>%
    as.matrix %>%
    as.data.frame %>%
    rownames_to_column('gene_id') %>%
    select(gene_id, everything()) %>%
    right_join(x=df_genes, by='gene_id') %>%
    write_tsv(snakemake@output$n_cells)
