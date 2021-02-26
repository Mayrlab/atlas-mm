configfile: "config.yaml"

import os

# make sure the tmp directory exists
os.makedirs(config["tmp_dir"], exist_ok=True)

# only use datasets we have in the config
datasets_selector = "(" + "|".join(config['sce']['txs'].keys()) + ")"


rule all:
    input:
        expand("data/sce/{dataset}.{level}.full_annot.Rds",
               dataset=['tmuris', 'brain', 'hspcs', 'merged'],
               level=['genes', 'txs']),
        "data/lui/merged_lui_expr_pointestimates.tsv.gz",
        "data/lui/merged_gene_expr_pointestimates.tsv.gz",
        "data/utrs/df_multiutrs.tsv"

rule merge_sces:
    input:
        tmuris=lambda wcs: config["sce"][wcs.level]["tmuris"],
        brain=lambda wcs: config["sce"][wcs.level]["brain"],
        hspcs=lambda wcs: config["sce"][wcs.level]["hspcs"]
    output:
        sce="data/sce/merged.{level}.raw.Rds"
    wildcard_constraints:
        level="(genes|txs)"
    conda:
        "envs/bioc_3_11.yaml"
    resources:
        mem_mb=12000
    script:
        "scripts/merge_sces.R"

rule compute_size_factors:
    input:
        sce="data/sce/merged.genes.raw.Rds"
    output:
        tsv="data/scran/merged.size_factors.tsv.gz"
    params:
        min_sf=0.95,
        max_sf=1.05,
        mrna_count=200000
    resources:
        mem_mb=16000
    conda:
        "envs/bioc_3_11.yaml"
    script:
        "scripts/compute_size_factors.R"

rule generate_utr_metadata:
    input:
        sce="data/sce/merged.txs.raw.Rds",
        ipa=config['utromeIPA']
    output:
        n_cells="data/utrs/n_cells_genes.tsv",
        utr="data/utrs/utr_metadata.tsv"
    params:
        min_cells=50
    resources:
        mem_mb=12000
    conda:
        "envs/bioc_3_11.yaml"
    script:
        "scripts/generate_utr_metadata.R"

rule generate_utr_length_metadata:
    input:
        utrome=config['utromeGTF'],
        gencode=config['gencodeValidEndsGTF'],
        utrs="data/utrs/utr_metadata.tsv",
        blacklist=config["utromeBlacklist"]
    output:
        txs="data/utrs/txs_utr_metadata_lengths.tsv",
        genes="data/utrs/genes_utr_metadata_lengths.tsv"
    resources:
        mem_mb=2000
    threads: 20
    conda:
        "envs/bioc_3_11.yaml"
    script:
        "scripts/generate_utr_length_metadata.R"

rule annotate_txs_sce:
    input:
        sce=lambda wcs: config['sce']['txs'][wcs.dataset],
        utrs="data/utrs/txs_utr_metadata_lengths.tsv",
        size_factors="data/scran/merged.size_factors.tsv.gz"
    output:
        sce="data/sce/{dataset}.txs.full_annot.Rds"
    wildcard_constraints:
        dataset=datasets_selector
    resources:
        mem_mb=4000
    conda:
        "envs/bioc_3_11.yaml"
    script:
        "scripts/annotate_txs_sce.R"

rule annotate_txs_sce_all:
    input:
        sce="data/sce/merged.txs.raw.Rds",
        utrs="data/utrs/txs_utr_metadata_lengths.tsv",
        size_factors="data/scran/merged.size_factors.tsv.gz"
    output:
        sce="data/sce/merged.txs.full_annot.Rds"
    resources:
        mem_mb=4000
    conda:
        "envs/bioc_3_11.yaml"
    script:
        "scripts/annotate_txs_sce_all.R"

rule annotate_genes_sce:
    input:
        sce=lambda wcs: config['sce']['genes'][wcs.dataset],
        utrs="data/utrs/genes_utr_metadata_lengths.tsv",
        size_factors="data/scran/merged.size_factors.tsv.gz"
    output:
        sce="data/sce/{dataset}.genes.full_annot.Rds"
    wildcard_constraints:
        dataset=datasets_selector
    resources:
        mem_mb=4000
    conda:
        "envs/bioc_3_11.yaml"
    script:
        "scripts/annotate_genes_sce.R"

rule annotate_sce_all:
    input:
        sce="data/sce/merged.genes.raw.Rds",
        utrs="data/utrs/genes_utr_metadata_lengths.tsv",
        size_factors="data/scran/merged.size_factors.tsv.gz"
    output:
        sce="data/sce/merged.genes.full_annot.Rds"
    resources:
        mem_mb=4000
    conda:
        "envs/bioc_3_11.yaml"
    script:
        "scripts/annotate_genes_sce_all.R"

rule generate_lui_table:
    input:
        sce="data/sce/merged.txs.full_annot.Rds",
        blacklist=config["utromeBlacklist"]
    params:
        min_cells=50
    output:
        lui="data/lui/merged_lui_expr_pointestimates.tsv.gz",
        n_cells="data/lui/merged_ncells_expr.tsv"
    conda:
        "envs/bioc_3_11.yaml"
    resources:
        mem_mb=12000
    script:
        "scripts/generate_lui_table.R"

rule generate_gene_table:
    input:
        sce="data/sce/merged.genes.full_annot.Rds",
        blacklist=config["utromeBlacklist"]
    output:
        tsv="data/lui/merged_gene_expr_pointestimates.tsv.gz"
    conda:
        "envs/bioc_3_11.yaml"
    resources:
        mem_mb=12000
    script:
        "scripts/generate_gene_table.R"

rule export_multiutrs:
    input:
        utrs="data/utrs/txs_utr_metadata_lengths.tsv",
        blacklist=config["utromeBlacklist"]
    output:
        tsv="data/utrs/df_multiutrs.tsv"
    conda:
        "envs/bioc_3_11.yaml"
    resources:
        mem_mb=4000
    script:
        "scripts/export_multiutrs.R"

rule export_multiutr_genes:
    input:
        multiutrs="data/utrs/df_multiutrs.tsv"
    output:
        ens="data/utrs/genes_multiutr_ensembl.txt",
        sym="data/utrs/genes_multiutr_symbols.txt",
        mgi="data/utrs/genes_multiutr_mgi.txt"
    conda:
        "envs/bioc_3_11.yaml"
    resources:
        mem_mb=4000
    script:
        "scripts/export_multiutr_genes.R"

rule genewalk_multiutr_genes:
    input:
        "data/utrs/genes_multiutr_mgi.txt"
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
