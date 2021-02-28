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
        input=list(utrs="data/utrs/utrome_genes_annotation.tsv"),
        output=list(ens="/fscratch/fanslerm/genes_multiutr_ensembl.txt",
                    sym="/fscratch/fanslerm/genes_multiutr_symbols.txt",
                    mgi="/fscratch/fanslerm/genes_multiutr_mgi.txt"),
        params=list())
}

################################################################################
## Load Data
################################################################################

## Load UTRs
df_genes <- read_tsv(snakemake@input$utrs) %>%
    filter(!is_blacklisted, atlas.utr_type == 'multi')

################################################################################
## Filter Overlapping Transcripts, Export
################################################################################

df_genes %>%
    pull(gene_id) %>%
    str_extract("^[^.]+") %>%
    unique() %>%
    write_lines(snakemake@output$ens)

df_genes %>%
    pull(gene_symbol) %>%
    unique() %>%
    write_lines(snakemake@output$sym)

mgis <- df_genes %>%
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
