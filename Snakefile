configfile: "config.yaml"

import pandas as pd
import os
from glob import glob

# make sure the tmp directory exists
os.makedirs(config["tmp_dir"], exist_ok=True)


rule all:
    input:
        expand("data/sce/merged.{level}.raw.rds", level=['genes', 'txs']),
        "data/scran/merged.size_factors.tsv.gz"


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
    threads: 8
    resources:
        mem=8
    conda:
        "envs/r36-scran.yaml"
    shell:
        """
        {input.script} {input.sce} {params.minSF} {params.maxSF} {params.mRNACount} {threads} {output}
        """
