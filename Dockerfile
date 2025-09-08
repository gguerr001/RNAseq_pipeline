# Lightweight base with micromamba for fast env solves
FROM mambaorg/micromamba:1.5.8

# Workdir
WORKDIR /workspace

# Copy repo contents
COPY . /workspace

# Install Snakemake in base (rule envs will be created on the fly)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      snakemake=7.* python=3.11 \
    && micromamba clean -a -y

# Make entrypoint executable
RUN chmod +x /workspace/entrypoint.sh

ENTRYPOINT ["/workspace/entrypoint.sh"]