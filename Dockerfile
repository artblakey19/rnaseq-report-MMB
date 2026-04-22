# Bulk RNA-seq downstream pipeline — container image.
#
# Multi-arch: linux/amd64 + linux/arm64 (Apple Silicon).
# Per-rule R/Python envs under workflow/envs/*.yaml are created by Snakemake
# on first run (--use-conda); they live in .snakemake/conda/ on the mounted
# project directory and are reused across runs.
#
# Build (single arch, local):
#   docker build -t bulk-rnaseq:latest .
#
# Build (multi-arch, requires buildx):
#   docker buildx create --use
#   docker buildx build --platform linux/amd64,linux/arm64 \
#       -t bulk-rnaseq:latest --load .
#
# Run (bind-mount project at /project, pass snakemake args as CMD):
#   docker run --rm -v "$PWD":/project bulk-rnaseq -c1
#
# The entrypoint (docker/entrypoint.sh) picks up the host UID/GID from
# /project's ownership and switches to it via gosu, and sets HOME=/tmp so
# conda/mamba have a writable cache. No `-u` or `-e HOME=/tmp` needed.

FROM mambaorg/micromamba:latest

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# System packages that Snakemake and conda env resolution expect,
# plus gosu for runtime UID switching in the entrypoint.
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        git \
        ca-certificates \
        gosu \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Base env: Snakemake + mamba frontend + Snakefile helpers (pandas, pyyaml).
# Quarto / R / Bioconductor packages are NOT installed here — each per-rule
# conda env under workflow/envs/*.yaml is resolved lazily at runtime.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        snakemake \
        conda \
        mamba \
        pandas \
        pyyaml \
 && micromamba clean --all --yes

# Bake the pipeline source into the image. Config and data stay on the host
# and are provided via bind-mount at /project.
WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow ./workflow
COPY --chown=$MAMBA_USER:$MAMBA_USER report   ./report

# Install entrypoint + dispatch. entrypoint.sh runs as root (the final USER
# below) so gosu can drop privileges to the host user derived from /project
# ownership. bulk-rnaseq-entry.sh then branches on the first arg: `init` →
# interactive config generator; anything else → snakemake pipeline.
USER root
COPY docker/entrypoint.sh          /usr/local/bin/entrypoint.sh
COPY docker/bulk-rnaseq-entry.sh   /usr/local/bin/bulk-rnaseq-entry.sh
RUN chmod +x /usr/local/bin/entrypoint.sh /usr/local/bin/bulk-rnaseq-entry.sh

# Expected runtime layout (bind-mount your project root here):
#   /project/config/{config.yaml, samples.tsv, contrasts.tsv}
#   /project/<counts>.tsv           (path comes from config.yaml)
#   /project/multiqc_data/
#   /project/.snakemake/conda/      (created automatically; reused on rerun)
WORKDIR /project

ENTRYPOINT ["/usr/local/bin/entrypoint.sh", "/usr/local/bin/bulk-rnaseq-entry.sh"]
CMD ["-c1"]
