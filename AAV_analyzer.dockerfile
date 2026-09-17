# AAVision container
# Base: micromamba (mamba-org) — a small, fast, fully free conda-compatible
# image. Replaces the deprecated continuumio/miniconda3 base. All dependencies
# (bioinformatics tools + Python libraries) are installed from env.yaml into the
# image's `base` environment.

FROM mambaorg/micromamba:2.8.1

LABEL maintainer="Icham Lahbib <icham.lahbib@kuleuven.be>"
LABEL description="AAV sequencing analysis pipeline (Nextflow workflow)"
LABEL version="1.0.0"
# Links the GHCR package to the source repository (shows under the repo's
# Packages and enables the repo to inherit the package's access settings).
LABEL org.opencontainers.image.source="https://github.com/Iclah/AAVision"

# ----------------------------------------------------------------------
# Build the conda environment from the pinned spec.
# The micromamba base image runs as a non-root user ($MAMBA_USER) and installs
# into the pre-existing `base` env. The COPY must be chown'd to that user.
# `micromamba install -n base -f env.yaml` resolves everything in one solve.
# ----------------------------------------------------------------------
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba clean --all --yes

# ----------------------------------------------------------------------
# Make the environment's tools available on PATH for NON-interactive shells.
# This is essential for Nextflow: it wraps each task in its own /bin/bash
# invocation and does NOT run the micromamba entrypoint activation, so relying
# on activation alone would leave minimap2/samtools/python off PATH. Setting
# PATH explicitly guarantees the tools resolve however the container is invoked.
# The base env lives at /opt/conda.
# ----------------------------------------------------------------------
ENV PATH=/opt/conda/bin:$PATH
ENV PYTHONPATH=/app

# ----------------------------------------------------------------------
# Application code. Copy as root and make it world-readable/executable so the
# code runs under any UID — Nextflow launches tasks with `-u $(id -u):$(id -g)`
# (see docker.runOptions in nextflow.config), so the runtime UID is the host
# user's, not $MAMBA_USER.
# ----------------------------------------------------------------------
USER root
COPY aav_analyze.py /app/aav_analyze.py
COPY aav_analyzer/ /app/aav_analyzer/
RUN chmod -R a+rX /app

# Convenience CLI shim. NOTE: we deliberately do NOT set an ENTRYPOINT.
# Nextflow (>=22.08) does not override a container ENTRYPOINT and wraps each task
# in a bash script; a hard ENTRYPOINT would be prepended to that wrapper and
# break execution. CMD keeps the image Nextflow-compatible.
RUN printf '#!/bin/bash\nexec python /app/aav_analyze.py "$@"\n' > /usr/local/bin/aav-analyze \
    && chmod +x /usr/local/bin/aav-analyze

# Data mount point / default working directory.
RUN mkdir -p /data && chmod a+rwx /data
WORKDIR /data

# Default command (shows help when run with no arguments).
CMD ["python", "/app/aav_analyze.py", "--help"]

# Build:    docker build -t aavision:dev -f AAV_analyzer.dockerfile .
# Nextflow: handled automatically (Nextflow supplies its own command wrapper).
