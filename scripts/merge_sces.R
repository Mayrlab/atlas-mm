#!/usr/bin/env Rscript

library(SingleCellExperiment)
library(S4Vectors)
library(tidyverse)
library(magrittr)

if (interactive()) {
    args <- c("/data/mayrc/data/tmuris/sce/utrome/TabulaMuris.genes.annot.rds",
              "/data/mayrc/data/ximerakis19/sce/utrome/brain.aging.genes.annot.rds",
              "/data/mayrc/data/mm-hspc-44k/sce/utrome/HSPCs.WT.genes.annot.rds",
              "/scratch/fanslerm/merged.genes.raw.rds")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 4) {
        stop("Incorrect number of arguments!\nUsage:\n> annotate_genes_sce.R <tmurisSCE> <brainSCE> <hspcsSCE> <outFile>\n")
    }
}
arg.tmurisRDS <- args[1]
arg.brainRDS <- args[2]
arg.hspcsRDS <- args[3]
arg.outFile <- args[4]

## Load SCEs
sce.tmuris <- readRDS(arg.tmurisRDS)
sce.brain  <- readRDS(arg.brainRDS)
sce.hspcs  <- readRDS(arg.hspcsRDS)

## Verify rownames are identical
stopifnot(all(rownames(sce.tmuris) == rownames(sce.brain)),
          all(rownames(sce.tmuris) == rownames(sce.hspcs)))

## Strip Row Data
rowData(sce.brain) <- NULL
rowData(sce.hspcs) <- NULL

########################################
## Adjust Column Data to Match
########################################

## Tabula Muris
colData(sce.tmuris) %<>%
    as.data.frame() %>%
    dplyr::rename(cell_type=cell_ontology_class, cluster=cluster.ids, sample=channel) %>%
    mutate(age='young') %>%
    select('cell', 'tissue', 'cell_type', 'cluster', 'sample', 'age',
           'phase', 'G1', 'S', 'G2M', 'G1.norm', 'S.norm', 'G2M.norm') %>%
    DataFrame()

## Brain
colData(sce.brain) %<>%
    as.data.frame() %>%
    select('cell', 'cell_type', 'sample_id', 'age_group',
           'phase', 'G1', 'S', 'G2M', 'G1.norm', 'S.norm', 'G2M.norm') %>%
    dplyr::rename(sample=sample_id, age=age_group) %>%
    mutate(tissue='Brain',
           cluster=as.integer(as.factor(cell_type))) %>%
    select('cell', 'tissue', 'cell_type', 'cluster', 'sample', 'age',
           'phase', 'G1', 'S', 'G2M', 'G1.norm', 'S.norm', 'G2M.norm') %>%
    DataFrame()

## HSPCs
colData(sce.hspcs) %<>%
    as.data.frame() %>%
    select('cell', 'clusters', 'sample', 'louvain_R',
           'phase', 'G1', 'S', 'G2M', 'G1.norm', 'S.norm', 'G2M.norm') %>%
    dplyr::rename(cell_type=clusters, cluster=louvain_R) %>%
    mutate(tissue='Bone Marrow (LSK,LK)',
           age='young') %>%
    select('cell', 'tissue', 'cell_type', 'cluster', 'sample', 'age',
           'phase', 'G1', 'S', 'G2M', 'G1.norm', 'S.norm', 'G2M.norm') %>%
    DataFrame()

sce <- cbind(sce.tmuris, sce.brain, sce.hspcs)
colnames(sce) <- sce$cell

saveRDS(sce, arg.outFile)
