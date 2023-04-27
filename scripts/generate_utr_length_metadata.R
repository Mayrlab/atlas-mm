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
    utrs.gr <- gr %>% filter(type == 'three_prime_UTR')
    ## three_prime_UTR no longer includes STOP
    bind_ranges(utrs.gr %>% filter_by_overlaps(stop.gr, maxgap=1),      # exon with STOP codon
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
        input=list(utrome="/data/mayrc/data/mca/gff/utrome.e30.t5.gc25.pas3.f0.9999.w500.gtf.gz",
                   gencode="/data/mayrc/data/mca/gff/gencode.vM25.mRNA_ends_found.gff3.gz",
                   utrs="data/utrs/utr_metadata.tsv",
                   blacklist="data/blacklist.utrome.genes.txt"),
        output=list(txs="/fscratch/fanslerm/txs_utr_metadata_lengths.tsv",
                    genes="/fscratch/fanslerm/genes_utr_metadata_lengths.tsv"),
        params=list(),
        threads=10L)
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
    dplyr::select(transcript_id, transcript_name, utr_name) %>%
    mutate(transcript_gencode=extractRefTx(transcript_id),
           offset=extractOffset(transcript_id))

df.utr.map <- left_join(df.utrs, df.utrome, by="transcript_id") %>%
    dplyr::select(transcript_id, transcript_gencode, utr_name, offset) 

gencode.pc_utrs.gr <- gencode.gr %>%
    filter(transcript_type == "protein_coding", (type == 'three_prime_UTR') | (type == 'stop_codon')) %>%
    mutate(exon_number=as.numeric(exon_number))

################################################################################
## Compute UTR lengths
################################################################################

df.utr.lengths <- mcmapply(
    function (transcript_id, transcript_gencode, utr_name, offset) {
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
             utr_name=utr_name,
             utr_length=ifelse(is.na(improper), NA, sum(widths) + offset),
             improper=improper)
    },
    transcript_id=df.utr.map$transcript_id,
    transcript_gencode=df.utr.map$transcript_gencode,
    utr_name=df.utr.map$utr_name,
    offset=df.utr.map$offset,
    mc.cores=snakemake@threads, USE.NAMES=FALSE
) %>% t() %>% as.data.frame()

df.utr.lengths.clean <- df.utr.lengths %>%
    mutate_at(c(1,2,3), c(as.character)) %>%
    mutate_at(4, as.integer) %>%
    mutate_at(5, as.logical)

df.utrs.final <- df.utrs %>%
    left_join(df.utr.lengths.clean, by='transcript_id')

################################################################################
## MANUAL CORRECTIONS
################################################################################

## Manually correcting the following UTRs
cat("The following UTRs are being manually corrected. See code for details.\n\n")
df.utrs.final %>%
    filter(improper) %>%
    dplyr::select(transcript_id, gene_name, utr_name, counts.total,
           utr.type.pct10.no_ipa, utr.pct.no_ipa,
           is_ipa, utr_length) %>%
    print(n=Inf)

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
    mutate(is_blacklisted=any(is.na(utr_length) & !is_ipa & utr.type.pct10.no_ipa == 'multi' &
                         ((utr.ncelltypes.pct10.no_ipa > 0) | (utr.pct.no_ipa >= 0.1)))) %>%
    ungroup() %>%
    mutate(is_blacklisted=is_blacklisted | (gene_name %in% genes_blacklist))


################################################################################
## Compute Gene Data
################################################################################

df.utr.genes <- df.utrs.final %>%
    filter(!is_ipa, (utr.pct.no_ipa >= 0.1) | (utr.ncelltypes.pct10.no_ipa > 0)) %>%
    group_by(gene_name, gene_id, gene.ncelltypes.cells50.no_ipa,
             utr.type.pct10.no_ipa,
             utr.count.pct10.no_ipa, is_blacklisted) %>%
    mutate(utr_type=ifelse(max(utr_length) == utr_length, "LU",
                    ifelse(min(utr_length) == utr_length, "SU", "MU"))) %>%
    pivot_wider(names_from=utr_type, values_from=utr_length) %>%
    summarise(SU=collapseSU(SU), MU=collapseMU(MU), LU=max(LU, na.rm=TRUE), .groups='drop') %>%
    dplyr::rename(ncelltypes50=gene.ncelltypes.cells50.no_ipa,
           utr_type=utr.type.pct10.no_ipa,
           utr_count=utr.count.pct10.no_ipa) %>%
    dplyr::select(gene_name, gene_id, utr_type, utr_count, SU, MU, LU,
                  ncelltypes50, is_blacklisted)

################################################################################
## Export Data 
################################################################################
write_tsv(df.utrs.final, snakemake@output$txs)
write_tsv(df.utr.genes, snakemake@output$genes)
