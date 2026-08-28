FROM mambaorg/micromamba:1.5.8

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
    && rm -rf /var/lib/apt/lists/*

RUN micromamba install -y -n base -c bioconda -c conda-forge \
        fastq-scan=0.4.1 \
        mash=2.3 \
        python=3.10 \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"