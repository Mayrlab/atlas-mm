#!/usr/bin/env Rscript

library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scran)

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(sce="data/sce/merged.genes.raw.Rds"),
        output=list(tsv="/fscratch/fanslerm/merged.size_factors.tsv.gz"),
        params=list(min_sf=0.95, max_sf=1.05, mrna_count=200000))
}

## Load SCE
sce <- readRDS(snakemake@input$sce)

## Only consider cell types with 21 cells (min for scran)
idx_cells <- colData(sce) %>%
    as.data.frame %>%
    add_count(tissue, cluster) %>%
    pull(n) %>%
    map_lgl(~ .x > 20) %>%
    which

sce <- sce[, idx_cells]

idx_min_genes <- which(rowMeans(counts(sce)) > 0.1)

## Compute size factors based on cell types
sizeFactors(sce) <- sce[idx_min_genes,] %>%
    calculateSumFactors(clusters=paste(sce$tissue, sce$cluster, sep='.'))
                      
df <- colData(sce) %>%
    as.data.frame %>%
    mutate(size_factor=sizeFactors(sce),
           cts=colSums(counts(sce)))

## Compute truncated capture efficiency
ce_centered <- df %>%
    filter(size_factor > snakemake@params$min_sf, size_factor < snakemake@params$max_sf) %>%
    pull(cts) %>%
    median() / snakemake@params$mrna_count

df %>%
    mutate(capture_efficiency=ce_centered*size_factor*100) %>%
    select('cell_id', 'size_factor', 'capture_efficiency') %>%
    write_tsv(snakemake@output$tsv)
