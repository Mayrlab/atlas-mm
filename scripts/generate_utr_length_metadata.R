#!/usr/bin/env Rscript

################################################################################
## Libraries
################################################################################

library(GenomicRanges)
library(GenomicFeatures)
library(tidyverse)
library(plyranges)
library(parallel)

################################################################################
## Methods
################################################################################

extractOffset <- function (tx) {
    ifelse(str_detect(tx, "-UTR"),
           tx %>% str_extract("(-|)[0-9]+$") %>% as.numeric(),
           0)
}

extractRefTx <- function (tx) {
    tx %>% str_extract("^ENSMUST[0-9]+\\.[0-9]+")
}

extractUTR3 <- function (gr) {
    stop.gr <- gr %>%
        filter(type == 'stop_codon') %>%
        filter(exon_number == min(exon_number))                         # use first when split
    utrs.gr <- gr %>% filter(type == 'UTR')
    bind_ranges(utrs.gr %>% filter_by_overlaps(stop.gr),                # exon with STOP codon
                utrs.gr %>% filter(exon_number > stop.gr$exon_number))  # downstream exons
}

collapseMU <- function (mus) {
    mus <- purrr::discard(mus, is.na)
    if (length(mus) == 0) {
        NA_character_
    } else if (length(mus) == 1) {
        as.character(mus)
    } else {
        paste(mus, collapse=";")
    }
}

collapseSU <- function (sus) {
    res <- suppressWarnings(min(sus, na.rm=TRUE))
    if (res == Inf) {
        NA
    } else {
        res
    }
}


################################################################################
## Mock Data
################################################################################

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots=c(input='list', output='list',
                                               params='list', threads='integer'))
    snakemake <- Snakemake(
        input=list(utrome="/data/mayrc/data/mca/gff/adult.utrome.e3.t200.f0.999.w500.gtf",
                   gencode="/data/mayrc/db/mm10/gencode.vM21.annotation.mRNA_ends_found.gtf.gz",
                   utrs="data/utrs/utr_metadata.tsv",
                   blacklist="data/blacklist.utrome.genes.txt"),
        output=list(txs="/fscratch/fanslerm/txs_utr_metadata_lengths.tsv",
                    genes="/fscratch/fanslerm/genes_utr_metadata_lengths.tsv"),
        params=list(),
        threads=3L)
}

################################################################################
## Load Data
################################################################################

utrome.gr  <- read_gff(snakemake@input$utrome)
gencode.gr <- read_gff(snakemake@input$gencode)
df.utrs <- read_tsv(snakemake@input$utrs)
genes_blacklist <- read_lines(snakemake@input$blacklist)

df.utrome <- utrome.gr %>%
    filter(type == 'transcript') %>%
    mcols() %>% as.data.frame() %>%
    dplyr::select(transcript_id, transcript_name) %>%
    mutate(transcript_gencode=extractRefTx(transcript_name),
           offset=extractOffset(transcript_name))

df.utr.map <- left_join(df.utrs, df.utrome, by="transcript_id") %>%
    dplyr::select(transcript_id, transcript_gencode, offset) 

gencode.pc_utrs.gr <- gencode.gr %>%
    filter(transcript_type == "protein_coding", (type == 'UTR') | (type == 'stop_codon')) %>%
    mutate(exon_number=as.numeric(exon_number))

################################################################################
## Compute UTR lengths
################################################################################

df.utr.lengths <- mcmapply(
    function (transcript_id, transcript_gencode, offset) {
        widths <- gencode.pc_utrs.gr %>%
            filter(.$transcript_id == transcript_gencode) %>%
            extractUTR3() %>%
            arrange(exon_number) %>%
            width()
        ## txs without annotated UTRs are marked 'NA'
        ## txs where the offset is upstream of the last exon is marked improper
        improper <- ifelse(length(widths) == 0, NA, tail(widths, n=1) + offset < 0)
        list(transcript_id=transcript_id,
             transcript_gencode=transcript_gencode,
             utr_length=ifelse(is.na(improper), NA, sum(widths) + offset),
             improper=improper)
    },
    transcript_id=df.utr.map$transcript_id,
    transcript_gencode=df.utr.map$transcript_gencode,
    offset=df.utr.map$offset,
    mc.cores=snakemake@threads, USE.NAMES=FALSE
) %>% t() %>% as.data.frame()

df.utr.lengths.clean <- df.utr.lengths %>%
    mutate_at(c(1,2), c(as.character)) %>%
    mutate_at(3, as.integer) %>%
    mutate_at(4, as.logical)

df.utrs.final <- df.utrs %>%
    left_join(df.utr.lengths.clean, by='transcript_id')

## DEBUGGING
## saveRDS(df.utrs.final, "data/scratch/20201111-df-utrs-final.Rds")
## saveRDS(utrome.gr, "data/scratch/20201111-utrome-gr.Rds")
## saveRDS(gencode.pc_utrs.gr, "data/scratch/20201111-gencode-pc_utrs-gr.Rds")

## if (interactive()) {
##     df.utrs.final <- readRDS("data/scratch/20201111-df-utrs-final.Rds")
##     utrome.gr <- readRDS("data/scratch/20201111-utrome-gr.Rds")
##     gencode.pc_utrs.gr <- readRDS("data/scratch/20201111-gencode-pc_utrs-gr.Rds")
## }

################################################################################
## MANUAL CORRECTIONS
################################################################################

## Manually correcting the following UTRs
cat("The following UTRs are being manually corrected. See code for details.\n\n")
df.utrs.final %>%
    filter(improper) %>%
    select(transcript_id, gene_symbol, counts.total,
           utr.type.pct10.no_ipa, utr.pct.no_ipa,
           is_ipa, utr_length)

#' Eif4a2.1
#' Intersecting the UTR of the GENCODE transcript with our UTR yields a correct length
idx.Eif4a2.1 <- which(df.utrs.final$transcript_id == 'Eif4a2.1')
df.utrs.final[idx.Eif4a2.1, 'utr_length'] <- GenomicRanges::intersect(
    utrome.gr %>% filter(transcript_id == 'Eif4a2.1'),
    gencode.pc_utrs.gr %>%
    filter(transcript_id == df.utrs.final[[idx.Eif4a2.1, 'transcript_gencode']])) %>%
    width()

#' Epm2aip1.1
#' The correct UTR length can be counted using the length of the long UTR and adding the offset
idx.Epm2aip1.1 <- which(df.utrs.final$transcript_id == 'Epm2aip1.1')
df.utrs.final[idx.Epm2aip1.1, 'utr_length'] <- (df.utrs.final %>% filter(transcript_id == 'Epm2aip1.3'))$utr_length + df.utr.map[[idx.Epm2aip1.1, 'offset']]

#' Golph3.1
#' intersect with exons from 3' UTR to get corrected length
idx.Golph3.1 <- which(df.utrs.final$transcript_id == 'Golph3.1')
df.utrs.final[idx.Golph3.1, 'utr_length'] <- GenomicRanges::intersect(
    utrome.gr %>% filter(transcript_id == 'Golph3.1'),
    gencode.pc_utrs.gr %>%
    filter(transcript_id == df.utrs.final[[idx.Golph3.1, 'transcript_gencode']]) %>%
    extractUTR3()) %>%
    width()

#' Hnrnpa2b1.2
#' intersect with exons from 3' UTR; still includes splicing, so need to sum
idx.Hnrnpa2b1.2 <- which(df.utrs.final$transcript_id == 'Hnrnpa2b1.2')
df.utrs.final[idx.Hnrnpa2b1.2, 'utr_length'] <- GenomicRanges::intersect(
    utrome.gr %>% filter(transcript_id == 'Hnrnpa2b1.2'),
    gencode.pc_utrs.gr %>%
    filter(transcript_id == df.utrs.final[[idx.Hnrnpa2b1.2, 'transcript_gencode']]) %>%
    extractUTR3()) %>%
    width() %>% sum()

#' Orai1.1
#' intersect with exons from 3' UTR
idx.Orai1.1 <- which(df.utrs.final$transcript_id == 'Orai1.1')
df.utrs.final[idx.Orai1.1, 'utr_length'] <- GenomicRanges::intersect(
    utrome.gr %>% filter(transcript_id == 'Orai1.1'),
    gencode.pc_utrs.gr %>%
    filter(transcript_id == df.utrs.final[[idx.Orai1.1, 'transcript_gencode']]) %>%
    extractUTR3()) %>%
    width()

#' Pigb.4
#' Compute distance in genomic coordinates from STOP codon to 3' end
idx.Pigb.4 <- which(df.utrs.final$transcript_id == 'Pigb.4')
df.utrs.final[idx.Pigb.4, 'utr_length'] <- (
    gencode.pc_utrs.gr %>%
    filter(transcript_id == df.utrs.final[[idx.Pigb.4, 'transcript_gencode']]) %>%
    filter(type == 'stop_codon') %>% end()) -
    (utrome.gr %>% filter(transcript_id == 'Pigb.4', type == 'exon') %>% start())

#' Primpol.2
#' intersect with exons from 3' UTR
idx.Primpol.2 <- which(df.utrs.final$transcript_id == 'Primpol.2')
df.utrs.final[idx.Primpol.2, 'utr_length'] <- GenomicRanges::intersect(
    utrome.gr %>% filter(transcript_id == 'Primpol.2'),
    gencode.pc_utrs.gr %>%
    filter(transcript_id == df.utrs.final[[idx.Primpol.2, 'transcript_gencode']]) %>%
    extractUTR3()) %>%
    width()

#' Prss35.1
#' Compute distance in genomic coordinates from STOP codon to 3' end
idx.Prss35.1 <- which(df.utrs.final$transcript_id == 'Prss35.1')
df.utrs.final[idx.Prss35.1, 'utr_length'] <- (
    utrome.gr %>% filter(transcript_id == 'Prss35.1', type == 'exon') %>% end()) -
    (gencode.pc_utrs.gr %>%
     filter(transcript_id == df.utrs.final[[idx.Prss35.1, 'transcript_gencode']]) %>%
     filter(type == 'stop_codon') %>% start())

#' Rbm3.2
#'
idx.Rbm3.2 <- which(df.utrs.final$transcript_id == 'Rbm3.2')
df.utrs.final[idx.Rbm3.2, 'utr_length'] <- GenomicRanges::intersect(
    utrome.gr %>% filter(transcript_id == 'Rbm3.2'),
    gencode.pc_utrs.gr %>%
    filter(transcript_id == df.utrs.final[[idx.Rbm3.2, 'transcript_gencode']]) %>%
    extractUTR3()) %>%
    width()


################################################################################
## Blacklist Multi-UTR Genes
################################################################################

#' We have three reasons for blacklisting genes:
#'  - gene collisions, where isoforms from multiple genes cannot be distinguished
#'  - downstream A-rich regions lead oversensitivity to fluctuations in internal
#'    priming, which we know occurs at a batch-level
#'  - unannotated "stop_codon" or "three_prime_utr", leading to zero length UTRs
#'
#' The former two are external and are maintained in the blacklist file. The
#' latter are computed here. All are marked here in the annotation.

## Identify genes where isoforms with unidentified UTRs are part of a multi-UTR gene
df.utrs.final %<>%
    group_by(gene_id) %>%
    mutate(is_blacklisted=any(is.na(utr_length) && !is_ipa && utr.type.pct10.no_ipa == 'multi' &&
                         ((utr.ncelltypes.pct10.no_ipa > 0) || (utr.pct.no_ipa >= 0.1)))) %>%
    ungroup() %>%
    mutate(is_blacklisted=is_blacklisted | (gene_symbol %in% genes_blacklist))

################################################################################
## Compute Gene Data
################################################################################

df.utr.genes <- df.utrs.final %>%
    filter(!is_ipa, (utr.pct.no_ipa >= 0.1) | (utr.ncelltypes.pct10.no_ipa > 0)) %>%
    group_by(gene_symbol, gene_id, gene.ncelltypes.cells50.no_ipa,
             utr.type.pct10.no_ipa,
             utr.count.pct10.no_ipa, is_blacklisted) %>%
    mutate(utr_type=ifelse(max(utr_length) == utr_length, "LU",
                    ifelse(min(utr_length) == utr_length, "SU", "MU"))) %>%
    spread(key=utr_type, value=utr_length) %>%
    summarise(SU=collapseSU(SU), MU=collapseMU(MU), LU=max(LU, na.rm=TRUE)) %>%
    ungroup() %>%
    rename(ncelltypes50=gene.ncelltypes.cells50.no_ipa,
           utr_type=utr.type.pct10.no_ipa,
           utr_count=utr.count.pct10.no_ipa) %>%
    dplyr::select(gene_symbol, gene_id, utr_type, utr_count, SU, MU, LU,
                  ncelltypes50, is_blacklisted)

################################################################################
## Export Data 
################################################################################
write_tsv(df.utrs.final, snakemake@output$txs)
write_tsv(df.utr.genes, snakemake@output$genes)
