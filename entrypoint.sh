#!/usr/bin/env bash
set -euo pipefail
cd /workspace

# Allow overrides
: "${CORES:=$(nproc)}"
: "${CONFIG:=config/config.yaml}"

snakemake -j "$CORES" --use-conda --conda-frontend mamba --configfile "$CONFIG" "$@"