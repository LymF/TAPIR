<div align="center">

<img src="logo.png" alt="TAPIR logo" width="400"/>

**Transcriptome Assembly Pipeline for Identification of RNA viruses**

[![Python](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.1.0-green.svg)]()
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)]()

</div>

---

## Overview

**TAPIR** is an end-to-end, checkpoint-aware pipeline for the discovery and annotation of RNA viruses from paired-end metatranscriptomics data. Starting from raw FASTQ files, TAPIR integrates quality control, host decontamination, dual-strategy *de novo* assembly, cross-assembly dereplication, contig extension, cross-sample consolidation, and taxonomic identification into a single, reproducible workflow.

TAPIR is designed for use with short paired-end Illumina reads and has been tested on metatranscriptomic data from environmental and host-associated samples.

---

## Pipeline overview

Steps 1–8 run independently for each sample. Steps 9–10 run once across all samples.

| Step | Tool | Description |
|------|------|-------------|
| 1 | fastp | Adapter trimming, quality filtering, PE error correction |
| 2 | Bowtie2 | Host genome decontamination — retain unmapped read pairs |
| 3a | rnaSPAdes | RNA-aware *de novo* assembly |
| 3b | SPAdes `--rnaviral` | RNA virus-optimised *de novo* assembly |
| 4 | MEGAHIT | Meta-sensitive *de novo* assembly |
| 5 | MMseqs2 | Pool all assemblers + dereplicate at 95% ANI (per sample) |
| 6 | Bowtie2 | Map cleaned reads back to assembly |
| 7 | CoverM | Per-contig mean coverage estimation |
| 8 | COBRA | Overlap-based contig extension |
| 9 | MMseqs2 | Cross-sample consolidation — concatenate all samples and dereplicate at 95% ANI |
| 10 | ViralQuest | BLAST + HMM + optional LLM annotation (one run over all samples) |

---

## Requirements

### System
- Linux (x86_64)
- Python ≥ 3.11
- ≥ 64 GB RAM (128+ GB recommended for large datasets)
- ≥ 500 GB disk space (databases included)

### Software dependencies

| Tool | Version tested | Purpose |
|---|---|---|
| [fastp](https://github.com/OpenGene/fastp) | ≥ 0.23 | QC and adapter trimming |
| [Bowtie2](https://github.com/BenLangmead/bowtie2) | ≥ 2.5 | Host removal + read mapping |
| [SAMtools](https://github.com/samtools/samtools) | ≥ 1.18 | BAM processing |
| [SPAdes](https://github.com/ablab/spades) (rnaSPAdes) | ≥ 3.15 | RNA-aware assembly |
| [MEGAHIT](https://github.com/voutcn/megahit) | ≥ 1.2.9 | Complementary assembly |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | ≥ 13 | Assembly dereplication |
| [CoverM](https://github.com/wwood/CoverM) | ≥ 0.6 | Coverage estimation (optional, has fallback) |
| [COBRA](https://github.com/linxingchen/cobra) (`cobra-meta`) | ≥ 1.2.3 | Contig extension |
| [ViralQuest](https://github.com/gabrielvpina/viralquest) | ≥ 0.1 | Viral identification |
| [Biopython](https://biopython.org/) | ≥ 1.81 | FASTA utilities |


---

## Installation

Three installation methods are available. All result in a `tapir` command available in your terminal.

---

### Option A — conda (recommended)

> Installs TAPIR and all external tools in one step.
> *(bioconda submission pending — use the manual method below until the package is available)*

```bash
# Once published to bioconda:
conda install -c bioconda -c conda-forge tapir-pipeline
tapir --help
```

**Manual conda install (available now):**

```bash
# 1. Clone the repository
git clone https://github.com/LymF/TAPIR.git
cd TAPIR

# 2. Create environment with all tools
mamba create -n tapir python=3.11 \
  -c bioconda -c conda-forge \
  fastp bowtie2 samtools \
  "spades>=3.15" megahit mmseqs2 coverm \
  --channel-priority flexible -y

conda activate tapir

# 3. Install Python dependencies and the tapir command
pip install cobra-meta viralquest biopython
pip install .

tapir --version
```

---

### Option B — Docker

> Fully self-contained — no environment setup required.
> Image available at: `ghcr.io/lymf/tapir:latest`

```bash
# Pull and run
docker pull ghcr.io/lymf/tapir:latest
docker run --rm -v /your/data:/data ghcr.io/lymf/tapir:latest \
    -i /data/reads -o /data/results \
    --host-genome /data/host.fa \
    -t 16 --ram 64 --email your@email.edu
```

**On HPC/shared servers without Docker root access — use Singularity/Apptainer:**

```bash
# Pull Docker image as a Singularity image file
singularity pull tapir.sif docker://ghcr.io/lymf/tapir:latest

# Run
singularity run tapir.sif \
    -i /data/reads -o /data/results \
    --host-genome /data/host.fa \
    -t 16 --ram 64 --email your@email.edu
```

> **Note:** Docker does not resolve AVX2 incompatibility — if the host CPU lacks AVX2,
> see [Tools on servers without AVX2](#tools-on-servers-without-avx2) below.

---

### Option C — pip

> Installs the `tapir` command and all Python dependencies. External bioinformatics tools (fastp, bowtie2, etc.) must be installed separately via conda.

```bash
pip install tapir-pipeline
tapir --help
```

> **Note:** If using a conda environment with Python ≠ 3.11, install biopython via conda first to avoid compilation errors:
> ```bash
> conda install -c bioconda biopython -y
> pip install tapir-pipeline
> ```

```bash
# Or install from source:
git clone https://github.com/LymF/TAPIR.git && cd TAPIR && pip install .
```

---

### Verify installation

```bash
tapir --version
# TAPIR 1.1.0
```

## Database setup

### RefSeq Viral (ViralQuest reference — ~219 MB)

```bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
gunzip viral.1.protein.faa.gz
diamond makedb --in viral.1.protein.faa --db viralDB.dmnd
```

### NCBI nr — DIAMOND format (~346 GB)

```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
diamond makedb --in nr --db nr.dmnd --threads 32
```

### HMM models

```bash
mkdir hmms && cd hmms

wget -O EggNOG-4.5.hmm.xz       https://zenodo.org/records/18715455/files/EggNOG-4.5.hmm.xz?download=1
wget -O U-RVDBv29.0-prot.hmm.xz https://zenodo.org/records/18715455/files/U-RVDBv29.0-prot.hmm.xz?download=1
wget -O Vfam-228.hmm.xz          https://zenodo.org/records/18715455/files/Vfam-228.hmm.xz?download=1
wget -O Pfam-A.hmm.xz            https://zenodo.org/records/18715455/files/Pfam-A.hmm.xz?download=1

unxz -v *.xz
```

---

## Usage

### Input format

Place paired-end FASTQ files in the input directory. Default naming convention:

```
/data/reads/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
└── sample2_R2.fastq.gz
```

Custom suffixes can be specified with `--r1-suffix` / `--r2-suffix`.

### Minimal run

```bash
tapir \
    -i /data/reads \
    -o /results \
    --host-genome /refs/host_genome.fa \
    -t 32 --ram 128 \
    --email your@email.edu
```

### Full run with all databases and LLM annotation

```bash
tapir \
    -i /data/reads \
    -o /results \
    --host-genome /refs/host_genome.fa \
    -t 64 --ram 256 \
    --email your@email.edu \
    --nr-dmnd    /dbs/nr.dmnd \
    --viral-dmnd /dbs/viralDB.dmnd \
    --rvdb-hmm   /dbs/hmms/U-RVDBv29.0-prot.hmm \
    --eggnog-hmm /dbs/hmms/eggNOG.hmm \
    --vfam-hmm   /dbs/hmms/Vfam228.hmm \
    --pfam-hmm   /dbs/hmms/Pfam-A.hmm \
    --llm-type google \
    --llm-model gemini-2.0-flash \
    --llm-api-key $GEMINI_KEY
```

### Skip host removal (pre-cleaned reads)

```bash
tapir -i /data/reads -o /results \
    --skip-host-removal \
    -t 32 --ram 128 --email your@email.edu
```

### Resume an interrupted run

TAPIR writes `.done_*` checkpoint files after each step. Re-run the same command to resume from the last successful step — no flags needed.

### Skip specific steps

```bash
tapir ... --skip-steps fastp host
# Available: fastp host rnaspades rnaviral megahit merge mapping coverage cobra cross_sample viralquest
```

### Local LLM via Ollama

```bash
tapir ... \
    --llm-type ollama \
    --llm-model qwen3:8b
# No API key required. Minimum recommended model: qwen3:4b
```

---

## Parameters reference

### Required

| Parameter | Description |
|---|---|
| `-i / --input-dir` | Directory containing paired FASTQ files |
| `-o / --output-dir` | Output directory |
| `--email` | Email address for NCBI online BLASTn |

### Resources

| Parameter | Default | Description |
|---|---|---|
| `-t / --threads` | `8` | CPU threads |
| `--ram` | `64` | Maximum RAM in GB |

### Host removal

| Parameter | Default | Description |
|---|---|---|
| `--host-genome` | — | Host reference genome FASTA |
| `--skip-host-removal` | `False` | Skip host decontamination |

### Assembly

| Parameter | Default | Description |
|---|---|---|
| `--mink` | `21` | Minimum k-mer size |
| `--maxk` | `141` | Maximum k-mer size (also sets COBRA expected overlap) |
| `--min-contig-len` | `500` | Minimum contig length after assembly |

### COBRA

| Parameter | Default | Description |
|---|---|---|
| `--cobra-query` | auto | Custom query FASTA; auto-selected if omitted |
| `--cobra-min-len` | `2000` | Minimum length for auto query selection |
| `--cobra-assembler` | `megahit` | Assembler hint for overlap calculation |

### Cross-sample consolidation (step 9)

| Parameter | Default | Description |
|---|---|---|
| `--cross-sample-id` | `0.95` | Min nucleotide identity for cross-sample MMseqs2 clustering |
| `--cross-sample-cov` | `0.95` | Min coverage of shorter sequence for cross-sample clustering |

### Databases (all optional but recommended)

| Parameter | Description |
|---|---|
| `--nr-dmnd` | DIAMOND nr database |
| `--viral-dmnd` | RefSeq viral DIAMOND database |
| `--rvdb-hmm` | RVDB protein HMM |
| `--eggnog-hmm` | eggNOG viral HMM |
| `--vfam-hmm` | Vfam HMM |
| `--pfam-hmm` | Pfam-A HMM |
| `--blastn-local PATH` | Local BLASTn database (overrides online BLASTn) |
| `--blastn-db DB` | NCBI nucleotide database for online BLASTn (default: `nt`) |
| `--max-orfs N` | Max non-overlapping ORFs per sequence for ViralQuest (default: `6`) |
| `--cap3` | Enable CAP3 contig assembly within ViralQuest (disabled by default) |

### LLM annotation

| Parameter | Description |
|---|---|
| `--llm-type` | Provider: `ollama` \| `openai` \| `anthropic` \| `google` |
| `--llm-model` | Model name (e.g. `gemini-2.0-flash`, `qwen3:8b`) |
| `--llm-api-key` | API key (required for cloud providers) |

---

## Output structure

At the end of the run TAPIR produces two output areas:

- **Per-sample directories** — full intermediate outputs for each sample (steps 1–8).
- **`final_results/`** — key deliverables organised into subfolders.

### `final_results/` — key deliverables

```
results/
└── final_results/
    ├── fastp_reports/
    │   ├── sample1_fastp.html          ← per-sample QC report
    │   ├── sample2_fastp.html
    │   └── ...
    ├── all_samples_viral.fa            ← final viral sequences (all samples)
    ├── all_samples_viral-BLAST.tsv     ← BLAST hit table (tab-separated)
    ├── all_samples_bestSeqs.json       ← per-sequence annotation (JSON)
    └── all_samples_visualization.html  ← interactive annotation report
```

Sequence headers in `all_samples_viral.fa` carry full provenance:
`>{sample}|{assembler}__{original_contig_id}`

### Full output tree

```
results/
├── tapir.log                           ← full pipeline log
├── final_results/                      ← see above
├── host_index/                         ← shared Bowtie2 host index (built once)
├── sample1/
│   ├── 01_fastp/
│   ├── 02_host_removal/
│   ├── 03_rnaspades/
│   ├── 04_megahit/
│   ├── 05_mmseqs/                      ← per-sample dereplication
│   ├── 06_mapping/
│   ├── 07_coverage/
│   └── 08_cobra/
├── sample2/  ...
├── 09_cross_sample/
│   └── all_samples_consolidated.fa     ← cross-sample dereplicated input to ViralQuest
└── 10_viralquest/
    └── all_samples/
        ├── all_samples_consolidated.fa_viral.fa
        ├── all_samples_consolidated.fa_viral-BLAST.csv
        ├── all_samples_consolidated.fa_bestSeqs.json
        └── all_samples_consolidated.fa_visualization.html
```

---

## Hardware recommendations

| Dataset size | Reads | CPU | RAM |
|---|---|---|---|
| Small | < 50 M | 16 | 64 GB |
| Medium | 50–200 M | 32 | 128 GB |
| Large | > 200 M | 64+ | 256+ GB |

> rnaSPAdes is the most RAM-intensive step. Reduce `--ram` if memory is limiting; SPAdes will stay within the budget at some cost to assembly quality.

---

## Checkpoint system

TAPIR writes a hidden `.done_<step>` sentinel file inside each step's output directory after successful completion. On a re-run the pipeline detects these flags and skips completed steps automatically.

- **Resume** an interrupted run: re-run the same command.
- **Re-run a step**: delete its `.done_*` file (e.g. `rm results/sample1/05_mmseqs/.done_merge`).
- **Re-run everything**: delete the output directory.

---

## Citation

If you use TAPIR in your research, please cite this repository and the tools it depends on:

**TAPIR pipeline**
> [Pending publication]

**COBRA**
> Chen, L., Banfield, J.F. COBRA improves the completeness and contiguity of viral genomes assembled from metagenomes. *Nat Microbiol* (2024). https://doi.org/10.1038/s41564-023-01598-2

**ViralQuest**
> Rodrigues, G.V.P., Ferreira, L.Y.M. & Aguiar, E.R.G.R. ViralQuest: a user-friendly interactive pipeline for viral-sequences analysis and curation. BMC Bioinformatics 27, 64 (2026). https://doi.org/10.1186/s12859-026-06391-6 — see https://github.com/gabrielvpina/viralquest

**SPAdes / rnaSPAdes**
> Prjibelski A. et al. Using SPAdes de novo assembler. *Curr Protoc Bioinformatics* (2020). https://doi.org/10.1002/cpbi.102

**MEGAHIT**
> Li D. et al. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics* (2015). https://doi.org/10.1093/bioinformatics/btv033

**MMseqs2**
> Steinegger M., Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat Biotechnol* (2017). https://doi.org/10.1038/nbt.3988

---

## Contributing

Contributions are welcome. Please open an issue to discuss proposed changes before submitting a pull request.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-improvement`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature/my-improvement`)
5. Open a Pull Request

---

## License

TAPIR is released under the [MIT License](LICENSE).

---

## Contact

For bug reports and feature requests, please use the [GitHub Issues](https://github.com/LymF/TAPIR/issues) page.

For general questions, contact: `lucasmelobiomed@gmail.com`
