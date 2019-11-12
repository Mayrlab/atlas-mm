#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(tidyverse)
library(org.Mm.eg.db)

################################################################################
## Load Arguments
################################################################################

if (interactive()) {
    args <- c("data/utrs/df-multiutrs.tsv",
              "/scratch/fanslerm/genes-multiutr-ensembl.txt",
              "/scratch/fanslerm/genes-multiutr-symbols.txt",
              "/scratch/fanslerm/genes-multiutr-mgi.txt")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 4) {
        stop("Incorrect number of arguments!\nUsage:\n> export_multiutr_genes.R <multiutrTSV> <outEnsemble> <outSymbols> <outMGI>\n")
    }
}
arg.utrs       <- args[1]
arg.outEnsembl <- args[2]
arg.outSymbols <- args[3]
arg.outMGI     <- args[4]

################################################################################
## Load Data
################################################################################

## Load UTRs
df.utrs <- read_tsv(arg.utrs)

################################################################################
## Filter Overlapping Transcripts, Export
################################################################################

df.genes <- df.utrs %>%
    filter(!overlapping) %>%
    group_by(gene_id) %>%
    mutate(tx_nonoverlapping=n()) %>%
    ungroup() %>%
    filter(tx_nonoverlapping > 1)

df.genes %>%
    pull(gene_id) %>%
    str_extract("^[^.]+") %>%
    unique() %>%
    write_lines(arg.outEnsembl)

df.genes %>%
    pull(gene_symbol) %>%
    unique() %>%
    write_lines(arg.outSymbols)

mgis <- df.genes %>%
    pull(gene_id) %>%
    str_extract("^[^.]+") %>%
    unique() %>%
    mapIds(x=org.Mm.eg.db, column='MGI', keytype='ENSEMBL', multiVals='list')

names(mgis) <- NULL
mgis %>%
    unlist() %>%
    na.omit() %>%
    unique() %>%
    str_extract("[0-9]+$") %>%
    write_lines(arg.outMGI)
