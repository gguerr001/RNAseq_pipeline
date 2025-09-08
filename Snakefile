###############################################
# RNA-seq Snakemake pipeline
# Steps: Trimmomatic → HISAT2 → samtools → featureCounts → merge → DESeq2
###############################################

import csv
import os
import pandas as pd

configfile: "config/config.yaml"

SAMPLES_DF = pd.read_csv(config["samplesheet"])  # columns: sample,fq1,fq2,condition
SAMPLES = SAMPLES_DF["sample"].tolist()

# Load contrasts (two columns: case,control)
with open(config["contrasts"], newline="") as f:
    CONTRASTS = [(row["case"], row["control"]) for row in csv.DictReader(f)]

def contrast_outputs():
    outs = []
    for case, control in CONTRASTS:
        outs.append(f"results/deseq2/{case}_vs_{control}.tsv")
        outs.append(f"results/deseq2/{case}_vs_{control}.xlsx")
    return outs

# Index handling
BUILD_INDEX = bool(config.get("hisat2", {}).get("genome_fasta", ""))
if BUILD_INDEX:
    INDEX_PREFIX = config.get("hisat2", {}).get("index_prefix", "ref/index/genome")
    IDX_FILES = [f"{INDEX_PREFIX}.{i}.ht2" for i in range(1,9)]
else:
    INDEX_PREFIX = config.get("hisat2", {}).get("index_prefix")
    IDX_FILES = []

rule all:
    input:
        # final DESeq2 artifacts + merged counts
        "results/counts/combined_counts.tsv",
        "results/deseq2/normalized_counts.tsv",
        "results/deseq2/PCA.pdf",
        "results/deseq2/heatmap_top50.pdf",
        contrast_outputs(),
        # per-sample counts
        expand("results/counts/{sample}.counts.txt", sample=SAMPLES),
        # optional built index
        IDX_FILES

########################
# Optional: build HISAT2 index
########################
if BUILD_INDEX:
    rule hisat2_index:
        input:
            fasta=config["hisat2"]["genome_fasta"]
        output:
            IDX_FILES
        threads: config.get("threads", {}).get("index", 8)
        conda:
            "envs/rna-map.yaml"
        shell:
            """
            hisat2-build -p {threads} {input.fasta} {INDEX_PREFIX}
            """

########################
# Trimming (Trimmomatic)
########################
rule trim:
    input:
        r1=lambda wildcards: SAMPLES_DF.set_index("sample").loc[wildcards.sample, "fq1"],
        r2=lambda wildcards: SAMPLES_DF.set_index("sample").loc[wildcards.sample, "fq2"],
    output:
        r1_paired="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r1_unpaired="results/trimmed/{sample}_R1.unpaired.fastq.gz",
        r2_paired="results/trimmed/{sample}_R2.trimmed.fastq.gz",
        r2_unpaired="results/trimmed/{sample}_R2.unpaired.fastq.gz",
    params:
        adapters=config.get("trimmomatic", {}).get("adapters", "auto"),
        extra="{params}"  # placeholder to keep formatters happy
    threads: config.get("threads", {}).get("trim", 8)
    conda:
        "envs/rna-map.yaml"
    shell:
        r"""
        set -euo pipefail
        ADAPTERS="{params.adapters}"
        if [ "$ADAPTERS" = "auto" ]; then
            # try to discover conda-installed adapters file
            ADAPTERS=$(dirname $(which trimmomatic))/../share/trimmomatic/adapters/TruSeq3-PE.fa
        fi
        trimmomatic PE -threads {threads} -phred33 \
            {input.r1} {input.r2} \
            {output.r1_paired} {output.r1_unpaired} \
            {output.r2_paired} {output.r2_unpaired} \
            ILLUMINACLIP:${{ADAPTERS}}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

########################
# Align (HISAT2) → sorted BAM + index
########################
rule align_sort_index:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz",
        idx=IDX_FILES if BUILD_INDEX else INDEX_PREFIX + ".1.ht2"
    output:
        bam="results/align/{sample}.sorted.bam",
        bai="results/align/{sample}.sorted.bam.bai"
    params:
        idx=INDEX_PREFIX
    threads: config.get("threads", {}).get("align", 12)
    conda:
        "envs/rna-map.yaml"
    shell:
        """
        set -euo pipefail
        hisat2 -p {threads} -x {params.idx} -1 {input.r1} -2 {input.r2} \
          | samtools view -b - \
          | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

########################
# featureCounts (gene-level)
########################
rule featurecounts:
    input:
        bam="results/align/{sample}.sorted.bam",
        bai="results/align/{sample}.sorted.bam.bai",
        gtf=config["gtf"]
    output:
        counts="results/counts/{sample}.counts.txt"
    threads: config.get("threads", {}).get("counts", 12)
    params:
        feature_type=config.get("featureCounts", {}).get("feature_type", "transcript"),
        attribute=config.get("featureCounts", {}).get("attribute", "gene_id"),
        extra=config.get("featureCounts", {}).get("extra", "-O")
    conda:
        "envs/rna-map.yaml"
    shell:
        """
        featureCounts -T {threads} \
          -t {params.feature_type} -g {params.attribute} \
          -a {input.gtf} -o {output.counts} {params.extra} {input.bam}
        """

########################
# Merge counts across samples
########################
rule merge_counts:
    input:
        counts=expand("results/counts/{sample}.counts.txt", sample=SAMPLES),
        samplesheet=config["samplesheet"]
    output:
        matrix="results/counts/combined_counts.tsv"
    conda:
        "envs/rna-map.yaml"
    script:
        "scripts/merge_counts.py"

########################
# DESeq2: stats + PCA + heatmap
########################
rule deseq2:
    input:
        counts="results/counts/combined_counts.tsv",
        samplesheet=config["samplesheet"],
        contrasts=config["contrasts"]
    output:
        pca="results/deseq2/PCA.pdf",
        heatmap="results/deseq2/heatmap_top50.pdf",
        norm="results/deseq2/normalized_counts.tsv",
        # dynamic per-contrast outputs
        tsvs=contrast_outputs()
    params:
        outdir="results/deseq2"
    conda:
        "envs/r-deseq2.yaml"
    shell:
        """
        Rscript scripts/deseq2.R \
          --counts {input.counts} \
          --samples {input.samplesheet} \
          --contrasts {input.contrasts} \
          --outdir {params.outdir}
        """