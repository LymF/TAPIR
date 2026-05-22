# Changelog

All notable changes to TAPIR will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [1.0.0] — 2026-05-22

### Added
- **Triple-assembler strategy**: SPAdes `--rnaviral` (step 3b) alongside rnaSPAdes and MEGAHIT, recovering more RNA virus contigs
- **Cross-sample consolidation** (step 9): per-sample COBRA outputs renamed `>SAMPLE|contigID`, merged and deduplicated with MMseqs2 `easy-linclust` (95% ANI, cov-mode 1) into a single multi-sample FASTA
- **Global ViralQuest** (step 10): single ViralQuest run over all samples combined; outputs in `10_viralquest/OUTPUT_all_samples/` and aggregated to `final_results/`
- `--cross-sample-id` and `--cross-sample-cov` CLI flags to tune cross-sample dereplication thresholds
- `rnaviral` and `cross_sample` added to `--skip-steps` choices
- TAPIR logo

### Fixed
- MMseqs2 clustering switched to `--cluster-mode 0` (set cover) to avoid off-by-one crash in parallel greedy mode
- Removed `--kmer-per-seq 80` from MMseqs2 calls (caused index overflow)
- `step_merge_assemblies` now skips missing or empty FASTA files gracefully
- Non-extended contigs no longer discarded: `merged_fa + cobra_fa` both forwarded to cross-sample consolidation

### Changed
- ViralQuest now runs once globally (post-consolidation) instead of per-sample
- Output directories renumbered: `09_cross_sample/`, `10_viralquest/`
- Installation documented with `--channel-priority flexible` for conda compatibility
- Source-compilation instructions added for servers without AVX2 (SPAdes, MEGAHIT, MMseqs2)

---

## [0.1.0] — 2025-01-01

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
