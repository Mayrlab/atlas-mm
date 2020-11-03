configfile: "config.yaml"

import os

# make sure the tmp directory exists
os.makedirs(config["tmp_dir"], exist_ok=True)

# only use datasets we have in the config
datasets_selector = "(" + "|".join(config['sce']['txs'].keys()) + ")"


rule all:
    input:
        expand("data/sce/{dataset}.{level}.full_annot.rds",
               dataset=['tmuris', 'brain', 'hspcs', 'merged'],
               level=['genes', 'txs']),
        "data/lui/merged-lui-expr-pointestimates.tsv",
        "data/lui/merged-gene-expr-pointestimates.tsv.gz",
        "data/utrs/df-multiutrs.tsv"

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
        mem_mb=12000
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
        mem_mb=16000
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
        mem_mb=12000
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
        txs="data/utrs/txs-utr-metadata-lengths.tsv",
        genes="data/utrs/genes-utr-metadata-lengths.tsv"
    resources:
        mem_mb=1000
    threads: 20
    conda:
        "envs/r36-plyranges.yaml"
    shell:
        """
        {input.script} {input.utrome} {input.gencode} {input.utrs} {threads} {output.txs} {output.genes}
        """

rule annotate_sce:
    input:
        sce=lambda wcs: config['sce'][wcs.level][wcs.dataset],
        utrs="data/utrs/{level}-utr-metadata-lengths.tsv",
        size_factors="data/scran/merged.size_factors.tsv.gz",
        script="scripts/annotate_{level}_sce.R"
    output:
        "data/sce/{dataset}.{level}.full_annot.rds"
    wildcard_constraints:
        dataset=datasets_selector
    resources:
        mem_mb=4000
    conda:
        "envs/r36-sce.yaml"
    shell:
        """
        {input.script} {input.sce} {input.utrs} {input.size_factors} {output}
        """

rule annotate_sce_all:
    input:
        sce=lambda wcs: ".".join(["data/sce/merged", wcs.level, "raw.rds"]),
        utrs="data/utrs/{level}-utr-metadata-lengths.tsv",
        size_factors="data/scran/merged.size_factors.tsv.gz",
        script="scripts/annotate_{level}_sce_all.R"
    output:
        "data/sce/merged.{level}.full_annot.rds"
    resources:
        mem_mb=4000
    conda:
        "envs/r36-sce.yaml"
    shell:
        """
        {input.script} {input.sce} {input.utrs} {input.size_factors} {output}
        """

rule generate_lui_table:
    input:
        script="scripts/generate_lui_table.R",
        sce="data/sce/merged.txs.full_annot.rds",
        blacklist=config["utromeBlacklist"],
        overlaps=config["utromeOverlaps"]
    params:
        minCells=50
    output:
        lui="data/lui/merged-lui-expr-pointestimates.tsv",
        ncells="data/lui/merged-ncells-expr.tsv"
    conda:
        "envs/r36-sce.yaml"
    resources:
        mem_mb=10000
    shell:
        """
        {input.script} {input.sce} {input.blacklist} {input.overlaps} {params.minCells} {output.lui} {output.ncells}
        """

rule generate_gene_table:
    input:
        script="scripts/generate_gene_table.R",
        sce="data/sce/merged.genes.full_annot.rds",
        blacklist=config["utromeBlacklist"]
    output:
        "data/lui/merged-gene-expr-pointestimates.tsv.gz"
    conda:
        "envs/r36-sce.yaml"
    resources:
        mem_mb=10000
    shell:
        """
        {input.script} {input.sce} {input.blacklist} {output}
        """

rule export_multiutrs:
    input:
        script="scripts/export_multiutrs.R",
        utrs="data/utrs/txs-utr-metadata-lengths.tsv",
        blacklist=config["utromeBlacklist"],
        overlaps=config["utromeOverlaps"]
    output:
        "data/utrs/df-multiutrs.tsv"
    conda:
        "envs/r36-sce.yaml"
    resources:
        mem_mb=4000
    shell:
        """
        {input.script} {input.utrs} {input.blacklist} {input.overlaps} {output}
        """

rule export_multiutr_genes:
    input:
        script="scripts/export_multiutr_genes.R",
        multiutrs="data/utrs/df-multiutrs.tsv"
    output:
        ens="data/utrs/genes-multiutr-ensembl.txt",
        sym="data/utrs/genes-multiutr-symbols.txt",
        mgi="data/utrs/genes-multiutr-mgi.txt"
    conda:
        "envs/r36-annot-mm.yaml"
    resources:
        mem_mb=4000
    shell:
        """
        {input.script} {input.multiutrs} {output}
        """

rule genewalk_multiutr_genes:
    input:
        "data/utrs/genes-multiutr-mgi.txt"
    output:
        "data/genewalk/multiutrs/genewalk_results.csv"
    params:
        base_dir=config['tmp_dir'] + '/genewalk',
        project="multiutr-atlas"
    conda:
        "envs/genewalk.yaml"
    resources:
        mem_mb=2000
    threads: 16
    shell:
        """
        genewalk --nproc {threads} --project {params.project} --genes {input} --id_type mgi_id --base_folder {params.base_dir}
        mkdir -p data/genewalk/{params.project}
        cp {params.base_dir}/{params.project}/* data/genewalk/multiutrs
        rm -rf {params.base_dir}/{params.project}
        """
