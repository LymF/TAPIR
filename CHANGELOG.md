# Changelog

All notable changes to TAPIR will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [0.1.0] — 2025-XX-XX

### Added
- Initial release of TAPIR
- Full end-to-end pipeline: fastp → Bowtie2 (host) → rnaSPAdes → MEGAHIT
  → MMseqs2 → Bowtie2 (mapping) → CoverM → COBRA → ViralQuest
- Automatic paired-end FASTQ discovery from input directory
- Checkpoint-based resume system (`.done_*` sentinel files)
- `--skip-steps` flag for bypassing individual pipeline stages
- CoverM (primary) and jgi_summarize_bam_contig_depths (fallback) for coverage
- ViralQuest LLM annotation support (Ollama, OpenAI, Anthropic, Google)
- ANSI colour logging to stdout + plain-text log file
- Full docstrings on all pipeline functions for publication-quality code
