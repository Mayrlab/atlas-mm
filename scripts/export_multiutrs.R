#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(tidyverse)


################################################################################
## Load Arguments
################################################################################

if (interactive()) {
    args <- c("data/utrs/txs-utr-metadata-lengths.tsv",
              "data/blacklist.utrome.txt",
              "data/utrs/adult.utrome.e3.t200.f0.999.w500.overlaps.tsv",
              "/scratch/fanslerm/df-multiutrs.tsv")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 4) {
        stop("Incorrect number of arguments!\nUsage:\n> export_multiutrs.R <utrMetadata> <blacklist> <overlaps> <outFile>\n")
    }
}
arg.utrs      <- args[1]
arg.blacklist <- args[2]
arg.overlaps  <- args[3]
arg.outFile   <- args[4]

################################################################################
## Load Data
################################################################################

## Load UTRs
df.utrs <- read_tsv(arg.utrs)

## Known blacklist
genes.blacklist <- read_lines(arg.blacklist)

## Load UTR intersections for blacklisting
df.intersections <- read_tsv(arg.overlaps, col_names=c("tx1", "tx2"), col_types="cc")


################################################################################
## Identify Multiutrs
################################################################################

df.utrs %>%
    filter(utr.count.pct10.no_ipa >= 2,
           is_ipa == FALSE,
           gene.ncelltypes.cells50.no_ipa >= 2,
           (utr.pct.no_ipa >= 0.1) | (utr.ncelltypes.pct10.no_ipa > 0)) %>%
    filter(!(gene_symbol %in% genes.blacklist)) %>%
    mutate(utr_position=as.integer(str_extract(transcript_name, "[0-9]+$"))) %>%
    group_by(gene_symbol) %>%
    mutate(is_LU=(utr_length == max(utr_length, na.rm=TRUE))) %>%
    mutate(overlapping={
        df.intersections %>%
            filter(tx1 == transcript_name[is_LU],
                   tx2 %in% transcript_name[!is_LU]) %>%
            nrow() > 0}) %>%
    ungroup() %>%
    arrange(gene_symbol, utr_position) %>%
    write_tsv(arg.outFile)
