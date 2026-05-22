FROM condaforge/miniforge3:latest

LABEL maintainer="Lucas Melo <lucasmelobiomed@gmail.com>"
LABEL description="TAPIR - Transcriptome Assembly Pipeline for Identification of RNA viruses"
LABEL version="1.0.0"

# Install build tools needed to compile Python C extensions (e.g. orfipy)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install all bioinformatics tools via mamba (faster solver, bioconda pre-configured)
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.11 \
    fastp \
    bowtie2 \
    samtools \
    spades \
    megahit \
    mmseqs2 \
    coverm \
    && mamba clean -afy

# Install Python dependencies and TAPIR via pip
RUN pip install --no-cache-dir tapir-pipeline

# Create working directory for user data
WORKDIR /data

ENTRYPOINT ["tapir"]
CMD ["--help"]
