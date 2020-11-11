#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(tidyverse)
library(org.Mm.eg.db)

################################################################################
## Mock Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(multiutrs="data/utrs/df_multiutrs.tsv"),
        output=list(ens="/fscratch/fansler/genes_multiutr_ensembl.txt",
                    sym="/fscratch/fansler/genes_multiutr_symbols.txt",
                    mgi="/fscratch/fansler/genes_multiutr_mgi.txt"),
        params=list())
}

################################################################################
## Load Data
################################################################################

## Load UTRs
df.utrs <- read_tsv(snakemake@input$multiutrs)

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
    write_lines(snakemake@output$ens)

df.genes %>%
    pull(gene_symbol) %>%
    unique() %>%
    write_lines(snakemake@output$sym)

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
    write_lines(snakemake@output$mgi)
