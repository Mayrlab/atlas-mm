#!/usr/bin/env Rscript

library(SingleCellExperiment)
library(S4Vectors)
library(tidyverse)
library(magrittr)

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(
            tmuris="/data/mayrc/data/tmuris/sce/tmuris.utrome.txs.Rds",
            brain="/data/mayrc/fanslerm/scutr-quant/data/sce/ximerakis19.utrome.txs.Rds",
            hspcs="/data/mayrc/data/mm-hspc-44k/sce/hspcs.utrome.txs.Rds",
            mescs="/data/mayrc/fanslerm/scutr-quant/data/sce/guo19.utrome.txs.Rds"),
        output=list(sce="/fscratch/fanslerm/merged.txs.raw.Rds"),
        params=list())
}

## Load SCEs
sce.tmuris <- readRDS(snakemake@input$tmuris)
sce.brain  <- readRDS(snakemake@input$brain)
sce.hspcs  <- readRDS(snakemake@input$hspcs)
sce.mescs  <- readRDS(snakemake@input$mescs)

## Verify rownames are identical
stopifnot(all(rownames(sce.tmuris) == rownames(sce.brain)),
          all(rownames(sce.tmuris) == rownames(sce.hspcs)),
          all(rownames(sce.tmuris) == rownames(sce.mescs)))

## Strip Row Data
rowData(sce.brain) <- NULL
rowData(sce.hspcs) <- NULL
rowData(sce.mescs) <- NULL

## Exclude Untyped Cells
sce.brain %<>% `[`(,!is.na(sce.brain$cell_type))

## Subset mESCs and MEFs
sce.mescs %<>% `[`(,sce.mescs$sample_id %in% c('mESC', 'MEF'))

########################################
## Adjust Column Data to Match
########################################

## Tabula Muris
colData(sce.tmuris) %<>%
    as.data.frame %>%
    dplyr::rename(cell_id=cell, cell_type=cell_ontology_class, cluster=cluster.ids, sample=channel) %>%
    mutate(age='young') %>%
    select('cell_id', 'tissue', 'cell_type', 'cluster', 'sample', 'age') %>%
    DataFrame()

## Brain
colData(sce.brain) %<>%
    as.data.frame %>%
    select('cell_id', 'cell_type', 'sample_id.x', 'age') %>%
    dplyr::rename(sample=sample_id.x) %>%
    mutate(tissue='Brain',
           cluster=as.integer(as.factor(cell_type))) %>%
    select('cell_id', 'tissue', 'cell_type', 'cluster', 'sample', 'age') %>%
    DataFrame()

## HSPCs
colData(sce.hspcs) %<>%
    as.data.frame %>%
    select('cell_id', 'clusters', 'sample', 'louvain_R') %>%
    dplyr::rename(cell_type=clusters, cluster=louvain_R) %>%
    mutate(tissue='Bone Marrow (LSK,LK)',
           age='young',
           cell_id=str_replace(cell_id, '^([ACGT]{16})_(.*)$', '\\2_\\1')) %>%
    select('cell_id', 'tissue', 'cell_type', 'cluster', 'sample', 'age') %>%
    DataFrame()

## mESCs
colData(sce.mescs) %<>%
    as.data.frame %>%
    select('cell_id', 'sample_id') %>%
    dplyr::rename(sample=sample_id) %>%
    mutate(tissue='Embryo',
           age='E13.5',
           cell_type=sample,
           cluster=as.integer(as.factor(cell_type))) %>%
    select('cell_id', 'tissue', 'cell_type', 'cluster', 'sample', 'age') %>%
    DataFrame()

sce <- cbind(sce.tmuris, sce.brain, sce.hspcs, sce.mescs)
colnames(sce) <- sce$cell_id

saveRDS(sce, snakemake@output$sce)
