#!/usr/bin/env Rscript

library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scran)

if (interactive()) {
    args <- c("data/sce/merged.genes.raw.rds",
              "0.95", "1.05", "200000",
              "/scratch/fanslerm/merged.capture.efficiency.centered.tsv.gz")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 5) {
        stop("Incorrect number of arguments!\nUsage:\n> collapse_sce_txs.R <sceRDS> <minSizeFactor> <maxSizeFactor> <medianmRNACount> <outFile>\n")
    }
}
arg.sceRDS  <- args[1]
arg.minSF   <- as.numeric(args[2])
arg.maxSF   <- as.numeric(args[3])
arg.mRNA    <- as.numeric(args[4])
arg.outFile <- args[5]

## Load SCE
sce <- readRDS(arg.sceRDS)

## Only consider cell types with 21 cells (min for scran)
idx.cells <- colData(sce) %>%
    as.data.frame() %>%
    group_by(tissue, cluster) %>%
    mutate(n=n()) %>%
    ungroup() %>%
    pull(n) %>%
    map_lgl(~ .x > 20) %>%
    which()

sce <- sce[, idx.cells]

idx.min.genes <- which(rowMeans(counts(sce)) > 0.1)

## Compute size factors based on cell types
sizeFactors(sce) <- counts(sce[idx.min.genes,]) %>%
    computeSumFactors(clusters=paste(sce$tissue, sce$cluster, sep='.'))
                      
df <- colData(sce) %>%
    as.data.frame() %>%
    mutate(size_factor=sizeFactors(sce),
           cts=colSums(counts(sce)))

## Compute truncated capture efficiency
ce.centered <- df %>%
    filter(size_factor > arg.minSF, size_factor < arg.maxSF) %>%
    pull(cts) %>%
    median() / arg.mRNA

df %>%
    mutate(capture_efficiency=ce.centered*size_factor*100) %>%
    select('cell', 'size_factor', 'capture_efficiency') %>%
    write_tsv(arg.outFile)
