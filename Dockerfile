FROM continuumio/miniconda3:latest

LABEL maintainer="Lucas Melo <lucasmelobiomed@gmail.com>"
LABEL description="TAPIR - Transcriptome Assembly Pipeline for Identification of RNA viruses"
LABEL version="1.0.0"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential wget curl git cmake \
    && rm -rf /var/lib/apt/lists/*

# Install all bioinformatics tools via conda/bioconda
RUN conda install -y -n base \
    -c conda-forge -c bioconda \
    --channel-priority flexible \
    python=3.11 \
    fastp \
    bowtie2 \
    samtools \
    spades \
    megahit \
    mmseqs2 \
    coverm \
    && conda clean -afy

# Install Python dependencies via pip
RUN pip install --no-cache-dir biopython cobra-meta viralquest

# Copy pipeline script and install as CLI tool
WORKDIR /opt/tapir
COPY tapir.py .
COPY pyproject.toml .
COPY README.md .
COPY LICENSE .
RUN pip install --no-cache-dir .

# Create working directory for user data
WORKDIR /data

ENTRYPOINT ["tapir"]
CMD ["--help"]
