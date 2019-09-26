configfile: "config.yaml"

import pandas as pd
import os
from glob import glob

# make sure the tmp directory exists
os.makedirs(config["tmp_dir"], exist_ok=True)


rule all:
    input:
        expand("data/sce/merged.{level}.raw.rds", level=['genes', 'txs']),
        "data/scran/merged.size_factors.tsv.gz",
        "data/utrs/gene-utr-metadata-lengths.tsv"


rule merge_sces:
    input:
        tmuris=lambda wcs: config["sce"][wcs.level]["tmuris"],
        brain=lambda wcs: config["sce"][wcs.level]["brain"],
        hspcs=lambda wcs: config["sce"][wcs.level]["hspcs"],
        script="scripts/merge_sces.R"
    output:
        "data/sce/merged.{level}.raw.rds"
    wildcard_constraints:
        level="(genes|txs)"
    conda:
        "envs/r36-sce.yaml"
    resources:
        mem=32
    shell:
        """
        {input.script} {input.tmuris} {input.brain} {input.hspcs} {output}
        """

rule compute_size_factors:
    input:
        sce="data/sce/merged.genes.raw.rds",
        script="scripts/compute_size_factors.R"
    output:
        "data/scran/merged.size_factors.tsv.gz"
    params:
        minSF=0.95,
        maxSF=1.05,
        mRNACount=200000
    resources:
        mem=32
    conda:
        "envs/r36-scran.yaml"
    shell:
        """
        {input.script} {input.sce} {params.minSF} {params.maxSF} {params.mRNACount} {output}
        """

rule generate_utr_metadata:
    input:
        sce="data/sce/merged.txs.raw.rds",
        ipa=config['utromeIPA'],
        script="scripts/generate_utr_metadata.R"
    output:
        ncells="data/utrs/ncells-genes.tsv",
        utr="data/utrs/utr-metadata.tsv"
    params:
        minCells=50
    resources:
        mem=32
    conda:
        "envs/r36-sce.yaml"
    shell:
        """
        {input.script} {input.sce} {input.ipa} {params.minCells} {output.ncells} {output.utr}
        """

rule generate_gene_utr_metadata:
    input:
        utrome=config['utromeGTF'],
        gencode=config['gencodeValidEndsGTF'],
        utrs="data/utrs/utr-metadata.tsv",
        script="scripts/generate_gene_utr_metadata.R"
    output:
        utr="data/utrs/utr-metadata-lengths.tsv",
        genes="data/utrs/gene-utr-metadata-lengths.tsv"
    resources:
        mem=4
    threads: 16
    conda:
        "envs/r36-plyranges.yaml"
    shell:
        """
        {input.script} {input.utrome} {input.gencode} {input.utrs} {threads} {output.utr} {output.genes}
        """



        
