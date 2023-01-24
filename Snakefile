"""Hevesi et al., 2023"""
import pandas as pd
from os import listdir, rename, getcwd
from os.path import join, basename, dirname, abspath
from pathlib import Path
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("7.20.0")

##### load config and sample sheets #####
configfile: "config.yaml"
samples = pd.read_table(config["samples"]).set_index("cellranger", drop=False)

container: "docker://continuumio/miniconda3:4.4.10"

##### target rules #####

shell.executable("/bin/bash")

rule all:
    input:
        expand("souporcell/{cellranger}/clusters.tsv",
                cellranger=samples["cellranger"]),
        expand("cellbender/{cellranger}/{cellranger}_output.h5",
                cellranger=samples["cellranger"])

##### load rules #####

CELLRANGER="source /home/etretiakov/src/cellranger-7.1.0/sourceme.bash && cellranger "

rule cellranger_count:
    input:
        sample=directory("fastq/{cellranger}"),
        idx=directory("mm10_optimized")
    output:
        summary="cellranger/{cellranger}/outs/web_summary.html",
        bam="cellranger/{cellranger}/outs/possorted_genome_bam.bam",
    params:
        sample="cellranger/{cellranger}"
    threads: 32
    shell:
        ("{CELLRANGER} count --include-introns true \
            --id={params.sample} \
            --transcriptome={input.idx} \
            --fastqs={input.sample} \
            --jobmode=local \
            --localcores={threads} ")

rule souporcell:
    input:
        bam="cellranger/{cellranger}/outs/possorted_genome_bam.bam",
        barcodes="cellranger/{cellranger}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        fasta="mm10_optimized/fasta/genome.fa"
    output:
        dir=directory("souporcell/{cellranger}"),
        clusters="souporcell/{cellranger}/clusters.tsv"
    params: 2
    container:
        "shub://wheaton5/souporcell"
    threads: 20
    shell:
        ("souporcell_pipeline.py \
            -i {input.bam} \
            -b {input.barcodes} \
            -f {input.fasta} \
            -t {threads} \
            -o {output.dir} \
            -k {params}")

rule cellbender:
    input:
        "cellranger/{cellranger}/outs/raw_feature_bc_matrix.h5"
    output:
        "cellbender/{cellranger}/{cellranger}_output.h5"
    params:
        ndroplets=lambda wildcards: samples["ndroplets"][wildcards.cellranger],
        ncells=lambda wildcards: samples["ncells"][wildcards.cellranger]
    container:
        "docker://etretiakov/cellbender:v0.0.1"
    threads: 20
    resources:
        nvidia_gpu=1
    shell:
        ("cellbender remove-background \
            --input {input} \
            --output {output} \
            --cuda \
            --expected-cells {params.ncells} \
            --total-droplets-included {params.ndroplets} \
            --fpr 0.01 \
            --epochs 150")
