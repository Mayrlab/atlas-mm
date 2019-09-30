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
## Load Arguments
################################################################################

if (interactive()) {
    args <- c("data/sce/merged.txs.full_annot.rds",
              "data/blacklist.utrome.txt",
              "data/utrs/adult.utrome.e3.t200.f0.999.w500.overlaps.tsv",
              "50",
              "/scratch/fanslerm/TM.lui.tsv",
              "/scratch/fanslerm/TM.ncells.tsv")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 6) {
        stop("Incorrect number of arguments!\nUsage:\n> generate_lui_table.R <sceRDS> <blacklist> <overlaps> <minCells> <outLUI> <outNCells>\n")
    }
}
arg.sceRDS    <- args[1]
arg.blacklist <- args[2]
arg.overlaps  <- args[3]
arg.minCells  <- as.integer(args[4])
arg.outLUI    <- args[5]
arg.outNCells <- args[6]

################################################################################
## Load Data
################################################################################

## Load UTR Counts
sce <- readRDS(arg.sceRDS)

## Known blacklist
genes.blacklist <- read_lines(arg.blacklist)

## Load UTR intersections for blacklisting
df.intersections <- read_tsv(arg.overlaps, col_names=c("tx1", "tx2"), col_types="cc")


################################################################################
## Identify Multiutrs
################################################################################

df.multiutrs <- rowData(sce) %>%
    as.data.frame() %>%
    filter(utr.count.no_ipa == 2,
           is_ipa == FALSE,
           gene.ncelltypes.cells50.no_ipa >= 2,
           (utr.pct.no_ipa >= 0.1) | (utr.ncelltypes.no_ipa > 0)) %>%
    filter(!(gene_symbol %in% genes.blacklist)) %>%
    group_by(gene_symbol) %>%
    mutate(overlapping={
        df.intersections %>%
            filter(tx1 == transcript_name[1],
                   tx2 == transcript_name[2]) %>%
            nrow() > 0}) %>%
    ungroup() %>%
    mutate(utr.position=as.integer(str_extract(transcript_name, "[0-9]+$"))) %>%
    arrange(gene_symbol, utr.position)

cat(sprintf("Removing **%d genes** due to overlapping transcripts.\n",
(df.multiutrs %>% filter(overlapping) %>% nrow())/2))

df.multiutrs %<>% filter(!overlapping)

idx.utr2 <- df.multiutrs %>% pull(transcript_name)

cat(sprintf("Considering **%d genes**.\n", length(idx.utr2)/2))

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

## M: (genes) x (utrs)
M.genes <- rowData(sce[idx.utr2,])$gene_symbol %>% fac2sparse()

## M: (distal UTRs) x (utrs)
M.LU <- sparseMatrix(i=1:dim(M.genes)[1], 2*(1:dim(M.genes)[1]), x=1, dims=dim(M.genes), dimnames=dimnames(M.genes))


################################################################################
## COMPUTE LUI AND GENE COUNTS 
################################################################################

## transcripts counts
cts.celltypes <- counts(sce[idx.utr2,]) %*% (M.celltypes / sce$size_factor.merged)

## mean gene counts
mu.gene.celltypes <- t(t(M.genes %*% cts.celltypes) / colSums(M.celltypes))

## number of cells in each cell type expressing a given gene
expr.gene.celltypes <- ((M.genes %*% counts(sce[idx.utr2,])) > 0) %*% M.celltypes

## Point Estimates for LUI for all (gene) x (cell type)
lui.point <- (M.LU %*% cts.celltypes) / (M.genes %*% cts.celltypes)

## Mask all gene-cell combinations
lui.point[which(expr.gene.celltypes < arg.minCells)] <- NA


################################################################################
## Wide table 
################################################################################

df.lui <- lui.point %>%
    as.matrix() %>% as.data.frame() %>%
    mutate(gene=rownames(.)) %>%
    select(gene, everything())

df.gene.expr <- mu.gene.celltypes %>%
    as.matrix() %>% as.data.frame() %>%
    mutate(gene=rownames(.)) %>%
    select(gene, everything())

df.full <- df.lui %>%
    left_join(df.gene.expr, by='gene', suffix=c(" LUI", " Mean Expression"))


################################################################################
## Export
################################################################################

write_tsv(df.full, arg.outLUI)

## data.frame(celltype=colnames(M.celltypes), cells=colSums(M.celltypes)) %>%
##     write_tsv("data/tmuris-cells-per-celltype.tsv")

expr.gene.celltypes %>%
    as.matrix() %>%
    as.data.frame() %>%
    mutate(gene=rownames(.)) %>%
    select(gene, everything()) %>%
    write_tsv(arg.outNCells)
