# RNA-seq Snakemake Pipeline (Docker + Conda)

End-to-end RNA-seq analysis using Trimmomatic → HISAT2 → samtools → featureCounts → DESeq2, fully reproducible with Snakemake, per-rule conda envs, and a Docker image.

## Quick start (Docker)

```bash
# 1) Clone this repository
# 2) Put your FASTQs under data/ and your GTF under ref/annotation.gtf
#    EITHER provide ref/genome.fa (pipeline will build the index),
#    OR point to an existing HISAT2 index via config/config.yaml.

# 3) Edit config/samples.csv and config/contrasts.csv

# 4) Build and run
docker build -t rnaseq-pipeline .
docker run --rm -v $(pwd):/workspace rnaseq-pipeline
# or limit cores
# docker run --rm -e CORES=8 -v $(pwd):/workspace rnaseq-pipeline
```

## Running locally (no Docker)

```bash
# Requires: conda/mamba installed
mamba create -n snakerna -c conda-forge -c bioconda snakemake=7 python=3.11 -y
conda activate snakerna
snakemake -j 8 --use-conda --conda-frontend mamba --configfile config/config.yaml
```

## Configuration

- **config/samples.csv** – sample sheet with columns:
  - `sample` (unique ID), `fq1`, `fq2` (paths to FASTQs), `condition` (group label)
- **config/contrasts.csv** – rows of `case,control` for DE comparisons
- **config/config.yaml** – paths for reference files, featureCounts options, threads

### Reference options
- **Option A (reproducible)**: Provide `ref/genome.fa` (or any path) and set `hisat2.genome_fasta`. The pipeline builds an index at `hisat2.index_prefix`.
- **Option B (existing index)**: Leave `genome_fasta: ""` and set `index_prefix` to your prebuilt HISAT2 prefix.

### Trimming adapters
`trimmomatic.adapters: auto` attempts to use conda’s bundled `TruSeq3-PE.fa`. Alternatively, point to `ref/adapters/TruSeq3-PE.fa`.

## Outputs
- `results/trimmed/` – trimmed reads (paired + unpaired)
- `results/align/` – sorted BAMs (`.bam` + `.bai`)
- `results/counts/` – per-sample featureCounts and `combined_counts.tsv`
- `results/deseq2/` –
  - `normalized_counts.tsv`
  - `PCA.pdf` – sample PCA from VST
  - `heatmap_top50.pdf` – top 50 variable genes
  - `<case>_vs_<control>.tsv/.xlsx` – DE tables with LFC shrinkage where available


## Reproducibility
- Each rule runs in its own frozen conda env (`envs/*.yaml`).
- The provided Dockerfile pins Snakemake and Python; rule tools are solved per env.

## Troubleshooting
- **Trimmomatic adapters**: if `auto` fails, set an explicit path to `TruSeq3-PE.fa`.
- **GTF attributes**: pipeline defaults to `gene_id`. Change in `config.yaml` if your GTF differs.

- **Index**: If you already have an index, set `genome_fasta: ""` and `index_prefix` accordingly.
