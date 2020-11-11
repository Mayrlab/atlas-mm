#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(tidyverse)

################################################################################
## Mock Arguments
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list', params='list'))
    snakemake <- Snakemake(
        input=list(utrs="data/utrs/txs_utr_metadata_lengths.tsv",
                   blacklist="data/blacklist.utrome.txt"),
        output=list(tsv="/fscratch/fanslerm/df_multiutrs.tsv"),
        params=list())
}

################################################################################
## Load Data
################################################################################

## Load UTRs
df.utrs <- read_tsv(snakemake@input$utrs)

## Known blacklist
genes.blacklist <- read_lines(snakemake@input$blacklist)

################################################################################
## Identify Multiutrs
################################################################################

df.utrs %>%
    filter(utr.count.pct10.no_ipa >= 2,
           is_ipa == FALSE,
           gene.ncelltypes.cells50.no_ipa >= 2,
           (utr.pct.no_ipa >= 0.1) | (utr.ncelltypes.pct10.no_ipa > 0)) %>%
    filter(!(gene_symbol %in% genes.blacklist)) %>%
    mutate(utr_position=as.integer(str_extract(transcript_id, "[0-9]+$"))) %>%
    group_by(gene_symbol) %>%
    mutate(is_LU=(utr_length == max(utr_length, na.rm=TRUE))) %>%
    ungroup() %>%
    arrange(gene_symbol, utr_position) %>%
    write_tsv(snakemake@output$tsv)
