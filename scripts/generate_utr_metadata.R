#!/usr/bin/env Rscript

## Libraries
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(S4Vectors)

## Methods
DataFrame.as.tibble <- function (df) {
    df %>%
        as.data.frame() %>%
        rownames_to_column()
}

tibble.as.DataFrame <- function (tbl) {
    tbl %>%
        column_to_rownames() %>%
        DataFrame()
}

## Credit: https://stackoverflow.com/a/12649993/570918
sdiv <- function(X, Y, names=dimnames(X)) {
    sX <- summary(X)
    sY <- summary(Y)
    sRes <- merge(sX, sY, by=c("i", "j"))
    sparseMatrix(i=sRes[,1], j=sRes[,2], x=sRes[,3]/sRes[,4],
                 dimnames=names)
}

dropInf <- function (M) {
    tmp <- M
    tmp[is.infinite(tmp)] <- 0
    drop0(tmp)
}

## Load Arguments
if (interactive()) {
    args <- c("data/sce/merged.txs.raw.rds",
              "data/ipa/adult.utrome.e3.t200.f0.999.w500.ipa.tsv",
              "50",
              "/scratch/fanslerm/ncells-genes.tsv",
              "/scratch/fanslerm/utr-metadata.tsv")
} else {
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) != 5) {
        stop("Incorrect number of arguments!\nUsage:\n> generate_utr_metadata.R <sceRDS> <ipaFile> <minCells> <outNCells> <outUTRs>\n")
    }
}
arg.sceRDS   <- args[1]
arg.ipas     <- args[2]
arg.minCells <- as.integer(args[3])
arg.outCells <- args[4]
arg.outUTRs  <- args[5]

## Load Data
sce <- readRDS(arg.sceRDS)

df.ipas <- read_tsv(arg.ipas, col_types='cc')

## Recompute UTR Info
rowData(sce) %<>%
    DataFrame.as.tibble() %>%
    mutate(counts.total=rowSums(counts(sce))) %>%
    group_by(gene_symbol) %>%
    mutate(counts.total.gene=sum(counts.total)) %>%
    mutate(utr.pct=counts.total/counts.total.gene,
           utr.rank=rank(desc(counts.total), ties.method="first")) %>%
    mutate(expressed.tx=(counts.total > 0), expressed.gene=(counts.total.gene > 0)) %>%
    mutate(utr.count.mca=n(), utr.count.merged=sum(expressed.tx)) %>%
    mutate(utr.type.mca=ifelse(utr.count.mca > 1, 'multi', 'single'),
           utr.type.merged=ifelse(utr.count.merged > 1, 'multi', 'single')) %>%
    mutate(utr.type.matches=(utr.type.mca == utr.type.merged)) %>%
    ungroup() %>%
    tibble.as.DataFrame()

################################################################################
## CELL-TYPE UTR ANALYSIS
################################################################################

## M: (genes) x (utrs)
gene_ids <- rowData(sce)$gene_id
names(gene_ids) <- rownames(sce)
M.genes <- gene_ids %>% fac2sparse()

rowData(sce)["is_ipa"] <- rownames(sce) %in% df.ipas$tx_name

## M: (genes) x (utrs)
M.genes.no_ipa <- drop0(t(t(M.genes) * as.integer(!rowData(sce)$is_ipa)))

## Percent Usage without IPA
cts.txs.total <- Matrix(rowSums(counts(sce)), ncol=1, sparse=TRUE)
cts.genes.no_ipa.total <- M.genes.no_ipa %*% cts.txs.total
pct.txs.no_ipa.total <- cts.txs.total %>%
    sdiv(t(M.genes.no_ipa) %*% cts.genes.no_ipa.total)

rowData(sce)["utr.pct.no_ipa"] <- as.numeric(pct.txs.no_ipa.total)
rowData(sce)["utr.count.total.pct05.no_ipa"] <- as.numeric(t(M.genes) %*% M.genes.no_ipa %*% (pct.txs.no_ipa.total > 0.05))
rowData(sce)["utr.count.total.pct10.no_ipa"] <- as.numeric(t(M.genes) %*% M.genes.no_ipa %*% (pct.txs.no_ipa.total > 0.10))


## M: (cells) x (tissue, cell_type, age)
M.groups <- colData(sce)[,c('tissue', 'cell_type', 'age')] %>%
    as.data.frame() %>%
    mutate(group=paste(str_replace_all(tissue, "_", " "),
                       cell_type, age, sep=', ')) %>%
    pull(group) %>%
    fac2sparse() %>%
    t()

## Percent Usage across cell types
ncells.genes.no_ipa.groups <- ((M.genes.no_ipa %*% counts(sce)) > 0) %*% M.groups
valid.genes.no_ipa.groups <- drop0(ncells.genes.no_ipa.groups >= arg.minCells)
rowData(sce)["gene.ncelltypes.cells50.no_ipa"] <- (t(M.genes) %*% rowSums(valid.genes.no_ipa.groups)) %>%
    as.numeric()

cts.txs.groups <- counts(sce) %*% M.groups
cts.genes.no_ipa.groups <- M.genes.no_ipa %*% cts.txs.groups
pct.txs.no_ipa.groups <- cts.txs.groups %>%
    sdiv(t(M.genes.no_ipa) %*% (cts.genes.no_ipa.groups * valid.genes.no_ipa.groups))


rowData(sce)["utr.ncelltypes.pct05.no_ipa"] <- rowSums(pct.txs.no_ipa.groups >= 0.05)
rowData(sce)["utr.ncelltypes.pct10.no_ipa"] <- rowSums(pct.txs.no_ipa.groups >= 0.1)
rowData(sce)["utr.count.celltypes.pct05.no_ipa"] <- t(M.genes) %*% M.genes.no_ipa %*% (rowSums(pct.txs.no_ipa.groups >= 0.05) > 0) %>% as.numeric()
rowData(sce)["utr.count.celltypes.pct10.no_ipa"] <- t(M.genes) %*% M.genes.no_ipa %*% (rowSums(pct.txs.no_ipa.groups >= 0.10) > 0) %>% as.numeric()


rowData(sce)["utr.count.pct05.no_ipa"] <- t(M.genes) %*% M.genes.no_ipa %*%
    ((rowData(sce)$utr.pct.no_ipa >= 0.05) | (rowData(sce)$utr.ncelltypes.pct05.no_ipa > 0)) %>%
    as.numeric()
rowData(sce)["utr.count.pct10.no_ipa"] <- t(M.genes) %*% M.genes.no_ipa %*%
    ((rowData(sce)$utr.pct.no_ipa >= 0.10) | (rowData(sce)$utr.ncelltypes.pct10.no_ipa > 0)) %>%
    as.numeric()

nzero.pct05 <- rowData(sce) %>%
    as.data.frame() %>%
    filter(!is.na(utr.count.pct05.no_ipa), utr.count.pct05.no_ipa == 0) %>%
    pull(gene_symbol) %>% unique() %>% length()

nzero.pct10 <- rowData(sce) %>%
    as.data.frame() %>%
    filter(!is.na(utr.count.pct10.no_ipa), utr.count.pct10.no_ipa == 0) %>%
    pull(gene_symbol) %>% unique() %>% length()

nnoexpr.pct05 <- rowData(sce) %>%
    as.data.frame() %>%
    filter(!is.na(utr.count.pct05.no_ipa), utr.count.pct05.no_ipa == 0, !expressed.gene) %>%
    pull(gene_symbol) %>% unique() %>% length()

nnoexpr.pct10 <- rowData(sce) %>%
    as.data.frame() %>%
    filter(!is.na(utr.count.pct10.no_ipa), utr.count.pct10.no_ipa == 0, !expressed.gene) %>%
    pull(gene_symbol) %>% unique() %>% length()

cat(sprintf("Found %i genes with no reads in non-IPA isoforms at 5%% threshold.\n", nzero.pct05))
cat(sprintf("Of those, %i (%0.1f%%) have no gene expression.\n\n",
            nnoexpr.pct05, 100*nnoexpr.pct05/nzero.pct05))

cat(sprintf("Found %i genes with no reads in non-IPA isoforms at 10%% threshold.\n", nzero.pct10))
cat(sprintf("Of those, %i (%0.1f%%) have no gene expression.\n\n",
            nnoexpr.pct10, 100*nnoexpr.pct10/nzero.pct10))

################################################################################
## EXPORT RESULTS
################################################################################

## Write Number of Cells Per Gene Per Cell Type Group
ncells.genes.no_ipa.groups %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column('gene_id') %>%
    write_tsv(arg.outCells)

rowData(sce) %>%
    as.data.frame() %>%
    mutate(utr.type.pct05.no_ipa=ifelse(utr.count.pct05.no_ipa > 1, "multi",
                                 ifelse(utr.count.pct05.no_ipa == 1, 'single',
                                        utr.type.mca)),
           utr.type.pct10.no_ipa=ifelse(utr.count.pct10.no_ipa > 1, "multi",
                                 ifelse(utr.count.pct10.no_ipa == 1, 'single',
                                        utr.type.mca))) %>%
    dplyr::rename(utr.pct.total=utr.pct, utr.rank.total=utr.rank) %>%
    select(transcript_name, gene_id, gene_symbol, chromosome,
           counts.total, counts.total.gene,
           utr.pct.total, utr.rank.total,
           expressed.tx, expressed.gene, is_ipa,
           utr.count.mca, utr.type.mca,
           utr.pct.no_ipa, utr.count.total.pct05.no_ipa, utr.count.total.pct10.no_ipa,
           utr.ncelltypes.pct05.no_ipa, utr.ncelltypes.pct10.no_ipa, gene.ncelltypes.cells50.no_ipa,
           utr.count.celltypes.pct05.no_ipa, utr.count.celltypes.pct10.no_ipa,
           utr.count.pct05.no_ipa, utr.count.pct10.no_ipa,
           utr.type.pct05.no_ipa, utr.type.pct10.no_ipa) %>%
      write_tsv(arg.outUTRs)
