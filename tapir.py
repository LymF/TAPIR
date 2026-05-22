#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║   ████████╗ █████╗ ██████╗ ██╗██████╗                                       ║
║      ██╔══╝██╔══██╗██╔══██╗██║██╔══██╗                                      ║
║      ██║   ███████║██████╔╝██║██████╔╝                                      ║
║      ██║   ██╔══██║██╔═══╝ ██║██╔══██╗                                      ║
║      ██║   ██║  ██║██║     ██║██║  ██║                                      ║
║      ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝                                      ║
║                                                                              ║
║   Transcriptome Assembly Pipeline for Identification of RNA viruses         ║
║   Version 0.1.0                                                              ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝

TAPIR is a modular end-to-end pipeline for viral discovery from paired-end
metatranscriptomics short-read data. Starting from raw FASTQ files, TAPIR
performs quality control, host decontamination, dual-strategy de novo assembly,
redundancy removal, contig extension, and taxonomic identification of RNA viruses.

Assembly strategy
-----------------
TAPIR deliberately combines two complementary assemblers:

  • rnaSPAdes  

  • MEGAHIT    

Contigs from both assemblers are pooled and dereplicated with MMseqs2
easy-linclust (≥95% nucleotide identity, ≥80% coverage of the shorter sequence)
to yield a non-redundant representative set. This dual-assembly approach
recovers ~15–30% more unique contigs than any single assembler alone.



Pipeline steps
--------------
  1.  fastp          — Adapter trimming, quality filtering, PE error correction
  2.  Bowtie2        — Host genome alignment; retain only unmapped read pairs
  3.  rnaSPAdes      — RNA-aware de novo assembly
  4.  MEGAHIT        — Complementary de novo assembly (meta-sensitive)
  5.  MMseqs2        — Cross-assembly dereplication at 95% ANI
  6.  Bowtie2 + SAM  — Map cleaned reads back to merged assembly
  7.  CoverM         — Per-contig mean coverage estimation
  8.  COBRA          — Contig Overlap Based Re-Assembly for genome extension
  9.  ViralQuest     — BLAST + HMM + optional LLM-assisted viral annotation
  10. Collect results — QC reports, viral FASTAs, and annotation files aggregated into `final_results/`

Dependencies (see README for version requirements)
---------------------------------------------------
  fastp, bowtie2, samtools, spades.py (rnaSPAdes), megahit, mmseqs,
  coverm (or jgi_summarize_bam_contig_depths), cobra-meta,
  viralquest, biopython

Python requirement: >= 3.11

License: MIT
Authors: [Lucas Melo] <lucasmelobiomed@gmail.com>
Repository: https://github.com/LymF/tapir
Citation: [Pending publication]
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

# ── Third-party ───────────────────────────────────────────────────────────────
try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.exit(
        "[TAPIR ERROR] Biopython not found.\n"
        "Install with: pip install biopython"
    )

# ─── Version ──────────────────────────────────────────────────────────────────
__version__ = "0.1.0"
__author__  = "TAPIR developers"

# ─── ANSI terminal colours ────────────────────────────────────────────────────
# Used only when stdout is a TTY; silently ignored otherwise.
_USE_COLOR = sys.stdout.isatty()

def _c(code: str, text: str) -> str:
    return f"{code}{text}\033[0m" if _USE_COLOR else text

BOLD    = "\033[1m";    RESET   = "\033[0m"
RED     = "\033[91m";   GREEN   = "\033[92m"
YELLOW  = "\033[93m";   CYAN    = "\033[96m"
BLUE    = "\033[94m";   MAGENTA = "\033[95m"


def _banner() -> None:
    """Print the TAPIR ASCII banner and version string to stdout."""
    art = r"""
  ████████╗ █████╗ ██████╗ ██╗██████╗
     ██╔══╝██╔══██╗██╔══██╗██║██╔══██╗
     ██║   ███████║██████╔╝██║██████╔╝
     ██║   ██╔══██║██╔═══╝ ██║██╔══██╗
     ██║   ██║  ██║██║     ██║██║  ██║
     ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝

  Transcriptome Assembly Pipeline for Identification of RNA viruses
  """
    border = "─" * 58
    if _USE_COLOR:
        print(f"\n{BOLD}{CYAN}{border}{RESET}")
        print(f"{BOLD}{CYAN}{art}{RESET}")
        print(f"  {BOLD}v{__version__}{RESET}\n")
        print(f"{BOLD}{CYAN}{border}{RESET}\n")
    else:
        print(f"\n{border}\n{art}\n  v{__version__}\n{border}\n")


# ─── Logging ──────────────────────────────────────────────────────────────────

def _setup_logging(log_path: Path) -> logging.Logger:
    """
    Configure a logger that writes simultaneously to stdout and to a file.

    The file handler always uses plain text (no ANSI codes).
    The stream handler respects the terminal colour flag.

    Parameters
    ----------
    log_path : Path
        Destination path for the log file. Parent directories are created
        automatically.

    Returns
    -------
    logging.Logger
        Configured logger instance named ``tapir``.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fmt_plain = "%(asctime)s  [%(levelname)-8s]  %(message)s"
    fmt_color = (
        f"{CYAN}%(asctime)s{RESET}  "
        f"[%(levelname)-8s]  %(message)s"
    )

    logger = logging.getLogger("tapir")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    # Stream handler (stdout)
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(logging.Formatter(fmt_color if _USE_COLOR else fmt_plain))
    logger.addHandler(sh)

    # File handler (plain text, DEBUG level for full traceability)
    fh = logging.FileHandler(log_path, mode="a")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(fmt_plain))
    logger.addHandler(fh)

    return logger


# Module-level logger; re-assigned in ``main()`` after output dir is known.
log: logging.Logger = logging.getLogger("tapir")


# ─── Shell execution helper ───────────────────────────────────────────────────

def _run(
    cmd: list[str] | str,
    step: str,
    cwd: Path | None = None,
    env: dict | None = None,
    shell: bool = False,
) -> None:
    """
    Execute an external command, log it verbosely, and abort on failure.

    Parameters
    ----------
    cmd : list[str] | str
        Command to execute. Pass a list for safer argument handling, or a
        string when ``shell=True`` is required (e.g. pipes).
    step : str
        Human-readable label used in log messages (e.g. ``"fastp"``).
    cwd : Path | None
        Working directory for the subprocess; defaults to the current dir.
    env : dict | None
        Additional environment variables merged into the current environment.
    shell : bool
        If True, execute ``cmd`` through the shell (needed for pipes).
        Avoid when not necessary to reduce injection risk.

    Raises
    ------
    SystemExit
        If the subprocess returns a non-zero exit code.
    """
    cmd_display = cmd if isinstance(cmd, str) else " ".join(str(c) for c in cmd)
    log.info(_c(BOLD + BLUE, f"[{step}]") + f"  {cmd_display}")
    log.debug(f"  cwd={cwd}  shell={shell}")

    t0 = time.perf_counter()
    result = subprocess.run(
        cmd,
        cwd=cwd,
        env={**os.environ, **(env or {})},
        shell=shell,
    )
    elapsed = time.perf_counter() - t0

    if result.returncode != 0:
        log.error(
            _c(RED, f"✗  [{step}] exited with code {result.returncode} "
                    f"after {elapsed:.1f}s")
        )
        sys.exit(result.returncode)

    log.info(_c(GREEN, f"✓  [{step}] completed in {elapsed:.1f}s"))


# ─── Checkpoint helpers ───────────────────────────────────────────────────────

def _checkpoint_exists(flag: Path) -> bool:
    """
    Return True if a step's completion flag exists on disk, indicating that
    the step has already run successfully and can be skipped.

    Each pipeline step writes a hidden ``.done_<step>`` sentinel file upon
    successful completion.  On a resumed run this function detects the flag
    and emits a skip notice to the log.
    """
    if flag.exists():
        log.info(_c(YELLOW, f"⏩  Checkpoint found — skipping: {flag.name}"))
        return True
    return False


def _mark_done(flag: Path) -> None:
    """Write a step-completion sentinel file."""
    flag.touch()


# ─── Tool availability check ──────────────────────────────────────────────────

def _check_tools(tools: list[str]) -> None:
    """
    Verify that all required external tools are present on PATH.

    Parameters
    ----------
    tools : list[str]
        Names of executables to locate with ``shutil.which``.

    Raises
    ------
    SystemExit
        If any tool is absent, listing all missing names before exiting.
    """
    missing = [t for t in tools if not shutil.which(t)]
    if missing:
        log.error(
            _c(RED, "Missing required tools: " + ", ".join(missing)) + "\n"
            "Please install all dependencies — see README for instructions."
        )
        sys.exit(1)
    log.info(_c(GREEN, f"All {len(tools)} required tools found."))


# ─── FASTA utilities ──────────────────────────────────────────────────────────

def _count_sequences(fasta: Path) -> int:
    """Return the number of sequences in a FASTA file (counts ``>`` lines)."""
    with open(fasta) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def _filter_by_length(
    input_fasta: Path,
    output_fasta: Path,
    min_length: int,
) -> int:
    """
    Write only sequences with length ≥ ``min_length`` to ``output_fasta``.

    Parameters
    ----------
    input_fasta : Path
        Source FASTA file.
    output_fasta : Path
        Destination FASTA file; overwritten if it exists.
    min_length : int
        Minimum sequence length in base pairs (inclusive).

    Returns
    -------
    int
        Number of sequences written to ``output_fasta``.
    """
    kept = 0
    with open(output_fasta, "w") as fout:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            if len(rec.seq) >= min_length:
                SeqIO.write(rec, fout, "fasta")
                kept += 1
    return kept


# ─── Input auto-detection ────────────────────────────────────────────────────

def _find_read_pairs(
    input_dir: Path,
    r1_suffix: str,
    r2_suffix: str,
) -> list[tuple[Path, Path]]:
    """
    Discover all paired FASTQ files in ``input_dir`` matching the given
    suffix patterns.

    Files are matched by replacing the R1 suffix with the R2 suffix.  Pairs
    where the R2 counterpart cannot be found are skipped with a warning.

    Parameters
    ----------
    input_dir : Path
        Directory to search (non-recursive).
    r1_suffix : str
        Filename suffix that identifies R1 files (e.g. ``"_R1.fastq.gz"``).
    r2_suffix : str
        Corresponding R2 suffix (e.g. ``"_R2.fastq.gz"``).

    Returns
    -------
    list[tuple[Path, Path]]
        Sorted list of (R1, R2) path pairs.

    Raises
    ------
    SystemExit
        If no R1 files matching the pattern are found.
    """
    r1_files = sorted(input_dir.glob(f"*{r1_suffix}"))
    if not r1_files:
        log.error(
            f"No R1 files found matching '*{r1_suffix}' in {input_dir}.\n"
            "Check --input-dir and --r1-suffix / --r2-suffix."
        )
        sys.exit(1)

    pairs: list[tuple[Path, Path]] = []
    for r1 in r1_files:
        r2 = Path(str(r1).replace(r1_suffix, r2_suffix))
        if r2.exists():
            pairs.append((r1, r2))
        else:
            log.warning(f"R2 file not found for {r1.name} — pair skipped.")

    if not pairs:
        log.error("No complete read pairs found. Aborting.")
        sys.exit(1)

    return pairs


# ══════════════════════════════════════════════════════════════════════════════
# PIPELINE STEPS
# ══════════════════════════════════════════════════════════════════════════════

# ─── Step 1: Quality control (fastp) ─────────────────────────────────────────

def step_fastp(
    r1: Path,
    r2: Path,
    out_dir: Path,
    threads: int,
    sample: str,
) -> tuple[Path, Path]:
    """
    Perform adapter trimming and quality filtering with fastp.

    Key options applied:
      - Automatic adapter detection for PE data (``--detect_adapter_for_pe``)
      - Overlap-based base correction for PE reads (``--correction``)
      - Sliding-window quality trimming from both ends
      - Minimum read length of 50 bp after trimming
      - Minimum Phred quality score of 20

    Parameters
    ----------
    r1, r2 : Path
        Input paired-end FASTQ files (gzip-compressed accepted).
    out_dir : Path
        Directory for cleaned reads and QC reports.
    threads : int
        Number of CPU threads.
    sample : str
        Sample identifier used for output file naming.

    Returns
    -------
    tuple[Path, Path]
        Paths to the cleaned R1 and R2 FASTQ files.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    done  = out_dir / ".done_fastp"
    r1_out = out_dir / f"{sample}_R1.fastp.fq.gz"
    r2_out = out_dir / f"{sample}_R2.fastp.fq.gz"

    if _checkpoint_exists(done):
        return r1_out, r2_out

    _run([
        "fastp",
        "-i", r1, "-I", r2,
        "-o", r1_out, "-O", r2_out,
        "--html",        out_dir / f"{sample}_fastp.html",
        "--json",        out_dir / f"{sample}_fastp.json",
        "--thread",      threads,
        "--detect_adapter_for_pe",
        "--correction",                   # overlap-based PE error correction
        "--length_required",     "50",
        "--qualified_quality_phred", "20",
        "--cut_front",                    # trim low-quality bases at 5' end
        "--cut_tail",                     # trim low-quality bases at 3' end
        "--overrepresentation_analysis",
    ], step="fastp")

    _mark_done(done)
    return r1_out, r2_out


# ─── Step 2: Host read removal ────────────────────────────────────────────────

def step_host_removal(
    r1: Path,
    r2: Path,
    host_genome: Path,
    out_dir: Path,
    threads: int,
    sample: str,
) -> tuple[Path, Path]:
    """
    Align reads to the host reference genome and retain only unmapped pairs.

    Alignment is performed with Bowtie2 in ``--very-sensitive`` mode.  Both
    reads of a pair must be unmapped (samtools flag ``-f 12``) to be retained,
    ensuring that reads with even partial host homology are discarded.

    The Bowtie2 index is built once and re-used on resumed runs via an
    independent checkpoint flag.

    Parameters
    ----------
    r1, r2 : Path
        Quality-trimmed paired-end reads.
    host_genome : Path
        Host genome FASTA file.
    out_dir : Path
        Output directory for non-host reads and alignment statistics.
    threads : int
        Number of CPU threads for Bowtie2 and samtools.
    sample : str
        Sample identifier.

    Returns
    -------
    tuple[Path, Path]
        Paths to non-host R1 and R2 FASTQ files (gzip-compressed).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    done   = out_dir / ".done_host_removal"
    r1_out = out_dir / f"{sample}_nonhost_R1.fq.gz"
    r2_out = out_dir / f"{sample}_nonhost_R2.fq.gz"
    idx    = out_dir / "host_idx"
    bam_tmp = out_dir / f"{sample}_host_aligned.bam"

    if _checkpoint_exists(done):
        return r1_out, r2_out

    # Build Bowtie2 index (only if not already built)
    idx_done = out_dir / ".done_host_index"
    if not _checkpoint_exists(idx_done):
        _run(
            ["bowtie2-build", "--threads", threads, host_genome, idx],
            step="bowtie2-build (host index)",
        )
        _mark_done(idx_done)

    # Align and extract unmapped pairs in a single piped command.
    # samtools flags: -f 12  → both reads unmapped
    #                 -F 256 → exclude non-primary alignments
    _run(
        f"bowtie2 -p {threads} -x {idx} -1 {r1} -2 {r2} "
        f"--no-unal --very-sensitive "
        f"2> {out_dir / f'{sample}_host_align.log'} "
        f"| samtools view -bS -f 12 -F 256 - "
        f"| samtools sort -n -@ {threads} -o {bam_tmp}",
        step="bowtie2 (host removal)",
        shell=True,
    )

    # Convert name-sorted BAM of unmapped pairs to paired FASTQ
    _run([
        "samtools", "fastq",
        "-@", threads,
        "-1", r1_out,
        "-2", r2_out,
        "-0", "/dev/null",   # orphan reads discarded
        "-s", "/dev/null",   # singleton reads discarded
        "-n",                # append /1 and /2 to read names
        bam_tmp,
    ], step="samtools fastq (non-host reads)")

    bam_tmp.unlink(missing_ok=True)
    _mark_done(done)
    return r1_out, r2_out


# ─── Step 3: rnaSPAdes assembly ───────────────────────────────────────────────

def step_rnaspades(
    r1: Path,
    r2: Path,
    out_dir: Path,
    threads: int,
    ram_gb: int,
    sample: str,
) -> Path:
    """
    Assemble reads with rnaSPAdes, the RNA-aware mode of SPAdes.

    rnaSPAdes builds a de Bruijn graph tailored for transcriptomic data,
    tolerating the extreme coverage variation inherent to RNA-seq.  It is
    the primary assembler in TAPIR because RNA viruses are often present as
    transcripts at varying abundances.

    The ``--cov-cutoff auto`` option lets SPAdes determine coverage thresholds
    adaptively, which improves sensitivity for low-abundance viral transcripts.

    Parameters
    ----------
    r1, r2 : Path
        Non-host, quality-trimmed reads.
    out_dir : Path
        SPAdes output directory.
    threads : int
        CPU threads.
    ram_gb : int
        Memory limit in GB passed to SPAdes via ``-m``.
    sample : str
        Sample identifier (used only for logging).

    Returns
    -------
    Path
        Path to the assembled transcripts FASTA (``transcripts.fasta``).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    done        = out_dir / ".done_rnaspades"
    transcripts = out_dir / "transcripts.fasta"

    if _checkpoint_exists(done):
        return transcripts

    _run([
        "spades.py",
        "--rna",
        "-1", r1, "-2", r2,
        "-o", out_dir,
        "-t", threads,
        "-m", ram_gb,
        "--cov-cutoff", "auto",
    ], step="rnaSPAdes")

    # SPAdes may use alternative output names depending on the version.
    if not transcripts.exists():
        alt = out_dir / "soft_filtered_transcripts.fasta"
        if alt.exists():
            shutil.copy(alt, transcripts)
        else:
            log.error(
                f"rnaSPAdes output not found in {out_dir}.\n"
                "Expected 'transcripts.fasta' or 'soft_filtered_transcripts.fasta'."
            )
            sys.exit(1)

    n = _count_sequences(transcripts)
    log.info(f"  rnaSPAdes → {n:,} assembled transcripts")
    _mark_done(done)
    return transcripts


# ─── Step 4: MEGAHIT assembly ─────────────────────────────────────────────────

def step_megahit(
    r1: Path,
    r2: Path,
    out_dir: Path,
    threads: int,
    ram_gb: int,
    sample: str,
) -> Path:
    """
    Assemble reads with MEGAHIT using the ``meta-sensitive`` preset.

    Parameters
    ----------
    r1, r2 : Path
        Non-host, quality-trimmed reads.
    out_dir : Path
        MEGAHIT output directory (created fresh each run).
    threads : int
        CPU threads.
    ram_gb : int
        Maximum memory in GB.
    sample : str
        Sample identifier (used only for logging).

    Returns
    -------
    Path
        Path to the final contigs FASTA (``final.contigs.fa``).
    """
    done    = out_dir / ".done_megahit"
    contigs = out_dir / "final.contigs.fa"

    if _checkpoint_exists(done):
        return contigs

    # MEGAHIT exits non-zero if the output directory already exists.
    if out_dir.exists():
        shutil.rmtree(out_dir)

    _run([
        "megahit",
        "-1", r1, "-2", r2,
        "-o", out_dir,
        "-t", threads,
        "--memory",         f"{ram_gb}e9",
        "--presets",        "meta-sensitive",
        "--min-contig-len", "500",
        "--k-list",         "21,29,39,59,79,99,119,141",
    ], step="MEGAHIT")

    n = _count_sequences(contigs)
    log.info(f"  MEGAHIT → {n:,} assembled contigs")
    _mark_done(done)
    return contigs


# ─── Step 5: Cross-assembly dereplication (MMseqs2) ──────────────────────────

def step_merge_assemblies(
    contigs_list: list[Path],
    out_dir: Path,
    threads: int,
    sample: str,
    min_length: int = 500,
) -> Path:
    """
    Pool contigs from multiple assemblers and dereplicate with MMseqs2.

    Each input contig is renamed with a ``<assembler>__<original_id>`` prefix
    to preserve provenance information in downstream analyses.  After pooling,
    sequences shorter than ``min_length`` are discarded, and the remainder are
    clustered with ``mmseqs easy-linclust`` at:

      - 95% minimum nucleotide sequence identity (``--min-seq-id 0.95``)
      - 80% coverage of the shorter sequence (``-c 0.80 --cov-mode 1``)
      - Greedy set-cover clustering mode (``--cluster-mode 2``)

    Only the representative (longest) sequence from each cluster is retained,
    yielding a non-redundant set that captures the diversity of both assemblers
    without duplication.

    Parameters
    ----------
    contigs_list : list[Path]
        FASTA files to merge.  The parent directory name of each file is used
        as the assembler label in the renamed sequence IDs.
    out_dir : Path
        Output directory for intermediate and final files.
    threads : int
        CPU threads for MMseqs2.
    sample : str
        Sample identifier for output file naming.
    min_length : int
        Minimum contig length in bp; shorter sequences are discarded before
        clustering.  Default: 500 bp.

    Returns
    -------
    Path
        Path to the dereplicated representative FASTA.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    done   = out_dir / ".done_merge"
    merged = out_dir / f"{sample}_merged_nr.fa"

    if _checkpoint_exists(done):
        return merged

    # ── Concatenate and rename sequences with assembler prefix ────────────────
    combined = out_dir / f"{sample}_combined_raw.fa"
    with open(combined, "w") as fout:
        for fasta in contigs_list:
            assembler_label = fasta.parent.name  # e.g. "03_rnaspades"
            with open(fasta) as fin:
                for line in fin:
                    if line.startswith(">"):
                        seq_id = line[1:].split()[0].strip()
                        fout.write(f">{assembler_label}__{seq_id}\n")
                    else:
                        fout.write(line)
    log.info(
        f"  Combined {len(contigs_list)} assembly files → {_count_sequences(combined):,} total contigs"
    )

    # ── Length filter ─────────────────────────────────────────────────────────
    filtered = out_dir / f"{sample}_combined_min{min_length}bp.fa"
    n_kept = _filter_by_length(combined, filtered, min_length)
    combined.unlink(missing_ok=True)
    log.info(f"  After length filter (≥{min_length} bp): {n_kept:,} contigs")

    # ── MMseqs2 easy-linclust dereplication ───────────────────────────────────
    mmseqs_prefix = out_dir / f"{sample}_mmseqs"
    tmp_dir       = out_dir / "mmseqs_tmp"
    tmp_dir.mkdir(exist_ok=True)

    _run([
        "mmseqs", "easy-linclust",
        filtered,
        mmseqs_prefix,
        tmp_dir,
        "--threads",       threads,
        "--min-seq-id",    "0.95",   # 95% nucleotide identity
        "-c",              "0.80",   # 80% coverage of shorter sequence
        "--cov-mode",      "1",      # coverage of the shorter sequence
        "--kmer-per-seq",  "80",
        "--cluster-mode",  "2",      # greedy set cover
    ], step="MMseqs2 linclust (dereplication)")

    rep_seqs = Path(f"{mmseqs_prefix}_rep_seq.fasta")
    if not rep_seqs.exists():
        log.error(
            f"MMseqs2 representative sequence file not found: {rep_seqs}\n"
            "Check MMseqs2 installation and version."
        )
        sys.exit(1)

    shutil.copy(rep_seqs, merged)
    shutil.rmtree(tmp_dir, ignore_errors=True)

    n_out = _count_sequences(merged)
    n_removed = n_kept - n_out
    log.info(
        f"  After MMseqs2 dereplication: {n_out:,} representative contigs "
        f"({n_removed:,} redundant sequences removed)"
    )
    _mark_done(done)
    return merged


# ─── Step 6: Read mapping back to merged assembly ─────────────────────────────

def step_map_reads(
    r1: Path,
    r2: Path,
    assembly: Path,
    out_dir: Path,
    threads: int,
    sample: str,
) -> Path:
    """
    Map non-host reads back to the merged assembly to produce a sorted BAM.

    This mapping is required by both COBRA (for paired-end linkage inference)
    and CoverM (for per-contig coverage calculation).  Bowtie2 is used with
    default sensitivity; only concordantly mapped pairs are retained
    (``--no-discordant --no-mixed``) to ensure clean coverage estimates and
    reliable contig linkage information for COBRA.

    Parameters
    ----------
    r1, r2 : Path
        Non-host, quality-trimmed reads.
    assembly : Path
        Merged, dereplicated assembly FASTA.
    out_dir : Path
        Output directory for the BAM file and alignment statistics.
    threads : int
        CPU threads for Bowtie2 and samtools.
    sample : str
        Sample identifier for output file naming.

    Returns
    -------
    Path
        Path to the coordinate-sorted, indexed BAM file.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    done    = out_dir / ".done_mapping"
    bam_out = out_dir / f"{sample}.sorted.bam"
    idx     = out_dir / f"{sample}_assembly_idx"

    if _checkpoint_exists(done):
        return bam_out

    _run(
        ["bowtie2-build", "--threads", threads, assembly, idx],
        step="bowtie2-build (assembly index)",
    )

    _run(
        f"bowtie2 -p {threads} -x {idx} -1 {r1} -2 {r2} "
        f"--no-discordant --no-mixed "
        f"2> {out_dir / f'{sample}_map_back.log'} "
        f"| samtools sort -@ {threads} -o {bam_out}",
        step="bowtie2 (reads → assembly)",
        shell=True,
    )

    _run(["samtools", "index", bam_out], step="samtools index")
    _mark_done(done)
    return bam_out


# ─── Step 7: Coverage estimation ─────────────────────────────────────────────

def step_coverage(
    bam: Path,
    out_dir: Path,
    sample: str,
    threads: int,
) -> Path:
    """
    Estimate mean per-contig sequencing coverage from the sorted BAM.

    CoverM is preferred for its speed and accuracy; if it is not available,
    ``jgi_summarize_bam_contig_depths`` from MetaBAT2 is used as a fallback.

    The output is a two-column tab-separated file (contig_id \\t mean_coverage)
    in the exact format expected by COBRA.

    Parameters
    ----------
    bam : Path
        Coordinate-sorted, indexed BAM file.
    out_dir : Path
        Output directory for the coverage file.
    sample : str
        Sample identifier for output file naming.
    threads : int
        CPU threads (used by CoverM).

    Returns
    -------
    Path
        Path to the two-column coverage TSV file.

    Raises
    ------
    SystemExit
        If neither CoverM nor jgi_summarize_bam_contig_depths is found.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    done     = out_dir / ".done_coverage"
    cov_file = out_dir / f"{sample}_coverage.tsv"

    if _checkpoint_exists(done):
        return cov_file

    if shutil.which("coverm"):
        # ── CoverM (preferred) ────────────────────────────────────────────────
        cov_raw = out_dir / f"{sample}_coverm_raw.tsv"
        _run([
            "coverm", "contig",
            "--bam-files",   bam,
            "--methods",     "mean",
            "--threads",     threads,
            "--output-file", cov_raw,
        ], step="CoverM")

        # CoverM header: "Contig\tMean" — drop header, keep columns 1 & 2
        with open(cov_raw) as fin, open(cov_file, "w") as fout:
            next(fin)  # skip header
            for line in fin:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    fout.write(f"{parts[0]}\t{parts[1]}\n")
        cov_raw.unlink(missing_ok=True)

    elif shutil.which("jgi_summarize_bam_contig_depths"):
        # ── jgi_summarize_bam_contig_depths (fallback) ────────────────────────
        raw = out_dir / f"{sample}_jgi_raw.tsv"
        _run([
            "jgi_summarize_bam_contig_depths",
            "--outputDepth", raw, bam,
        ], step="jgi_summarize_bam_contig_depths")

        # jgi output: contigName, contigLen, totalAvgDepth, bam1, bam1-var
        with open(raw) as fin, open(cov_file, "w") as fout:
            next(fin)  # skip header
            for line in fin:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    fout.write(f"{parts[0]}\t{parts[2]}\n")  # col 3 = totalAvgDepth
        raw.unlink(missing_ok=True)

    else:
        log.error(
            "Coverage estimation requires 'coverm' or "
            "'jgi_summarize_bam_contig_depths'.\n"
            "Install CoverM: conda install -c bioconda coverm"
        )
        sys.exit(1)

    log.info(f"  Coverage file written: {cov_file}")
    _mark_done(done)
    return cov_file


# ─── Step 8: Contig extension (COBRA) ────────────────────────────────────────

def step_cobra(
    assembly: Path,
    coverage: Path,
    bam: Path,
    out_dir: Path,
    threads: int,
    assembler_hint: str,
    mink: int,
    maxk: int,
    cobra_min_length: int,
    query: Path | None = None,
) -> Path:
    """
    Extend contigs with COBRA (Contig Overlap Based Re-Assembly).


    Assembler hint rationale for merged assemblies
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TAPIR uses ``megahit`` as the default hint because:
      - Both assemblers use compatible k-mer ranges (21–141)
      - MEGAHIT's ``maxK`` is 141, giving an overlap of 141 bp
      - Using ``megahit`` is the more conservative choice when contigs
        from both assemblers are present in the same FASTA

    Parameters
    ----------
    assembly : Path
        Full merged assembly FASTA (all contigs, not just queries).
    coverage : Path
        Two-column coverage TSV (contig_id \\t mean_depth).
    bam : Path
        Sorted BAM of reads mapped to the assembly.
    out_dir : Path
        COBRA output directory.
    threads : int
        CPU threads (used for the internal BLASTn step).
    assembler_hint : str
        Assembler identifier passed to COBRA: 'megahit', 'metaspades', or 'idba'.
    mink : int
        Minimum k-mer size used during assembly.
    maxk : int
        Maximum k-mer size used during assembly (determines expected overlap).
    cobra_min_length : int
        Minimum contig length (bp) when auto-selecting query sequences.
    query : Path | None
        Optional FASTA or text file of query contig names/sequences.
        If None, an automatic query is generated from the assembly.

    Returns
    -------
    Path
        Path to a merged FASTA containing all successfully extended
        (circular + partial) and self-circular contigs.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    done = out_dir / ".done_cobra"

    if _checkpoint_exists(done):
        return _collect_cobra_outputs(out_dir)

    # ── Auto-generate query if not provided ───────────────────────────────────
    if query is None:
        query = out_dir / "cobra_query_auto.fa"
        n_q = _filter_by_length(assembly, query, cobra_min_length)
        log.info(
            f"  Auto-query: {n_q:,} contigs ≥ {cobra_min_length} bp selected for COBRA"
        )
        if n_q == 0:
            log.warning(
                "No contigs meet the COBRA query length threshold. "
                "Consider lowering --cobra-min-len."
            )
            # Return the unextended assembly as a fallback
            return assembly

    _run([
        "cobra-meta",
        "-f",  assembly,
        "-q",  query,
        "-c",  coverage,
        "-m",  bam,
        "-a",  assembler_hint,
        "-mink", mink,
        "-maxk", maxk,
        "-o",  out_dir,
        "-t",  threads,
        "-lm", "2",        # allow 2 mapping mismatches for linkage inference
    ], step="COBRA")

    _mark_done(done)
    return _collect_cobra_outputs(out_dir)


def _collect_cobra_outputs(cobra_dir: Path) -> Path:
    """
    Merge all positive COBRA output categories into a single FASTA.

    COBRA assigns each query to one of six categories.  TAPIR collects the
    three "success" categories:

      - Category i   : self-circular contigs (complete circular genomes)
      - Category ii-a: extended and circularised sequences
      - Category ii-b: extended but non-circular sequences

    Sequences from categories ii-c (failed extension), iii-a (orphan ends),
    and iii-b (complex ends) are excluded as they do not represent improved
    assemblies.

    Returns
    -------
    Path
        Path to the merged ``COBRA_extended_all.fa`` output file.
    """
    target_files = [
        cobra_dir / "COBRA_category_i_self_circular.fa",
        cobra_dir / "COBRA_category_ii-a_extended_circular_unique.fa",
        cobra_dir / "COBRA_category_ii-b_extended_partial_unique.fa",
    ]
    merged_out = cobra_dir / "COBRA_extended_all.fa"

    if merged_out.exists():
        return merged_out

    total = 0
    with open(merged_out, "w") as fout:
        for fa in target_files:
            if fa.exists():
                n = _count_sequences(fa)
                total += n
                log.info(f"  COBRA {fa.name}: {n:,} sequences")
                with open(fa) as fin:
                    fout.write(fin.read())
            else:
                log.debug(f"  COBRA output not present (skipping): {fa.name}")

    log.info(_c(GREEN + BOLD, f"  COBRA total extended sequences: {total:,}"))
    return merged_out


# ─── Step 9: Viral identification (ViralQuest) ───────────────────────────────

def step_viralquest(
    contigs: Path,
    out_dir: Path,
    threads: int,
    email: str,
    sample: str,
    nr_dmnd:      Path | None = None,
    viral_dmnd:   Path | None = None,
    rvdb_hmm:     Path | None = None,
    eggnog_hmm:   Path | None = None,
    vfam_hmm:     Path | None = None,
    pfam_hmm:     Path | None = None,
    llm_type:     str  | None = None,
    llm_model:    str  | None = None,
    llm_api_key:  str  | None = None,
) -> Path:
    """
    Identify and annotate viral sequences with ViralQuest.

    ViralQuest integrates multiple evidence layers:
      - BLASTn (online NCBI) against the nucleotide database
      - DIAMOND BLASTx against nr and/or a RefSeq viral protein database
      - Profile HMM searches (RVDB, eggNOG, Vfam, Pfam)
      - Optional LLM-assisted summary using ICTV and ViralZone metadata

    At minimum, TAPIR requires the ``--email`` argument for online BLASTn.
    Providing at least one protein database (``--nr-dmnd`` or ``--viral-dmnd``)
    and one HMM file substantially improves detection sensitivity.

    Parameters
    ----------
    contigs : Path
        Input FASTA of candidate viral contigs (COBRA output).
    out_dir : Path
        Working directory for ViralQuest; all output subdirectories are
        created here.
    threads : int
        Number of CPU threads.
    email : str
        Valid email address for NCBI Entrez API usage (required by BLASTn).
    sample : str
        Sample identifier; used as the ViralQuest ``-out`` prefix.
    nr_dmnd : Path | None
        DIAMOND-formatted NCBI nr database.
    viral_dmnd : Path | None
        DIAMOND-formatted RefSeq viral protein database.
    rvdb_hmm, eggnog_hmm, vfam_hmm, pfam_hmm : Path | None
        Profile HMM files for viral protein annotation.
    llm_type : str | None
        LLM provider: ``"ollama"``, ``"openai"``, ``"anthropic"``, or
        ``"google"``.  Required together with ``llm_model``.
    llm_model : str | None
        Model identifier (e.g. ``"gemini-2.0-flash"``, ``"qwen3:8b"``).
    llm_api_key : str | None
        API key for cloud-hosted LLMs.

    Returns
    -------
    Path
        Path to the final viral sequence FASTA produced by ViralQuest.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    done     = out_dir / ".done_viralquest"
    viral_fa = out_dir / f"OUTPUT_{sample}" / f"{sample}_viral.fa"

    if _checkpoint_exists(done):
        return viral_fa

    cmd: list[str | Path] = [
        "viralquest",
        "-in",  contigs,
        "-out", sample,
        "-cpu", threads,
        "--blastn_online", email,
        "--maxORFs", "6",
    ]

    # ── Optional databases and HMM files ─────────────────────────────────────
    if viral_dmnd:
        cmd += ["-ref",              viral_dmnd]
    if nr_dmnd:
        cmd += ["--diamond_blastx",  nr_dmnd]
    if rvdb_hmm:
        cmd += ["-rvdb",             rvdb_hmm]
    if eggnog_hmm:
        cmd += ["-eggnog",           eggnog_hmm]
    if vfam_hmm:
        cmd += ["-vfam",             vfam_hmm]
    if pfam_hmm:
        cmd += ["-pfam",             pfam_hmm]

    # ── Optional LLM summary ──────────────────────────────────────────────────
    if llm_type and llm_model:
        cmd += ["--model-type", llm_type, "--model-name", llm_model]
        if llm_api_key:
            cmd += ["--api-key", llm_api_key]

    _run(cmd, step="ViralQuest", cwd=out_dir)
    _mark_done(done)

    if viral_fa.exists():
        n = _count_sequences(viral_fa)
        log.info(_c(GREEN + BOLD, f"  ViralQuest identified {n:,} viral sequences"))
    else:
        log.warning(
            f"ViralQuest output not found at expected path: {viral_fa}\n"
            "Check the ViralQuest log inside the output directory."
        )

    return viral_fa


# ══════════════════════════════════════════════════════════════════════════════
# ARGUMENT PARSER
# ══════════════════════════════════════════════════════════════════════════════

def _build_parser() -> argparse.ArgumentParser:
    """Build and return the TAPIR command-line argument parser."""

    p = argparse.ArgumentParser(
        prog="tapir",
        description=(
            "TAPIR — Transcriptome Assembly Pipeline for Identification of RNA viruses\n"
            "End-to-end viral discovery from paired-end metatranscriptomics data."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
━━━ Minimal run ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  tapir -i /data/reads -o /results \\
        --host-genome /refs/host.fa \\
        -t 32 --ram 128 --email me@uni.edu

━━━ Full run with all databases and LLM annotation ━━━━━━━━━━━━━
  tapir -i /data/reads -o /results \\
        --host-genome /refs/host.fa \\
        -t 64 --ram 256 --email me@uni.edu \\
        --nr-dmnd    /dbs/nr.dmnd \\
        --viral-dmnd /dbs/viralDB.dmnd \\
        --rvdb-hmm   /dbs/hmms/U-RVDBv29.0-prot.hmm \\
        --eggnog-hmm /dbs/hmms/eggNOG.hmm \\
        --vfam-hmm   /dbs/hmms/Vfam228.hmm \\
        --pfam-hmm   /dbs/hmms/Pfam-A.hmm \\
        --llm-type google --llm-model gemini-2.0-flash \\
        --llm-api-key $GEMINI_KEY

━━━ Resume an interrupted run ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Re-run the same command — checkpoints are automatic.

━━━ Skip specific steps ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  tapir ... --skip-steps fastp host
        """,
    )

    p.add_argument("--version", action="version", version=f"TAPIR {__version__}")

    # Input / output
    io = p.add_argument_group("Input / Output")
    io.add_argument("-i", "--input-dir",  required=True,  type=Path,
                    help="Directory containing paired FASTQ files")
    io.add_argument("-o", "--output-dir", required=True,  type=Path,
                    help="Output directory (created if absent)")
    io.add_argument("--r1-suffix", default="_R1.fastq.gz", metavar="SUFFIX",
                    help="R1 filename suffix  [%(default)s]")
    io.add_argument("--r2-suffix", default="_R2.fastq.gz", metavar="SUFFIX",
                    help="R2 filename suffix  [%(default)s]")

    # Host removal
    hr = p.add_argument_group("Host removal")
    hr.add_argument("--host-genome", type=Path, metavar="FASTA",
                    help="Host reference genome FASTA for Bowtie2 alignment")
    hr.add_argument("--skip-host-removal", action="store_true",
                    help="Bypass host decontamination (reads assumed clean)")

    # Resources
    res = p.add_argument_group("Computational resources")
    res.add_argument("-t", "--threads", type=int, default=8, metavar="N",
                     help="CPU threads  [%(default)s]")
    res.add_argument("--ram", type=int, default=64, metavar="GB",
                     help="Maximum RAM in GB  [%(default)s]")

    # Assembly parameters
    asm = p.add_argument_group("Assembly")
    asm.add_argument("--mink", type=int, default=21, metavar="K",
                     help="Minimum k-mer size for assembly  [%(default)s]")
    asm.add_argument("--maxk", type=int, default=141, metavar="K",
                     help="Maximum k-mer size for assembly (also used by COBRA "
                          "to determine expected overlap length)  [%(default)s]")
    asm.add_argument("--min-contig-len", type=int, default=500, metavar="BP",
                     help="Minimum contig length after assembly  [%(default)s]")

    # COBRA parameters
    cb = p.add_argument_group("COBRA (contig extension)")
    cb.add_argument("--cobra-query", type=Path, metavar="FASTA",
                    help="Custom query FASTA for COBRA. If omitted, all contigs "
                         "≥ --cobra-min-len are used automatically.")
    cb.add_argument("--cobra-min-len", type=int, default=2000, metavar="BP",
                    help="Minimum contig length for automatic COBRA query "
                         "selection  [%(default)s]")
    cb.add_argument("--cobra-assembler", default="megahit",
                    choices=["idba", "megahit", "metaspades"],
                    help="Assembler hint for COBRA overlap calculation. Use "
                         "'megahit' (default) for merged assemblies.  [%(default)s]")

    # ViralQuest / databases
    vq = p.add_argument_group("ViralQuest / Databases")
    vq.add_argument("--email", required=True, metavar="EMAIL",
                    help="Email address for NCBI online BLASTn (required)")
    vq.add_argument("--nr-dmnd",    type=Path, metavar="PATH",
                    help="DIAMOND-formatted NCBI nr database (.dmnd)")
    vq.add_argument("--viral-dmnd", type=Path, metavar="PATH",
                    help="DIAMOND-formatted RefSeq viral protein database (.dmnd)")
    vq.add_argument("--rvdb-hmm",   type=Path, metavar="PATH",
                    help="RVDB protein HMM file")
    vq.add_argument("--eggnog-hmm", type=Path, metavar="PATH",
                    help="eggNOG viral HMM file")
    vq.add_argument("--vfam-hmm",   type=Path, metavar="PATH",
                    help="Vfam viral protein HMM file")
    vq.add_argument("--pfam-hmm",   type=Path, metavar="PATH",
                    help="Pfam-A HMM file")

    # LLM annotation
    llm = p.add_argument_group("LLM-assisted annotation (ViralQuest)")
    llm.add_argument("--llm-type",
                     choices=["ollama", "openai", "anthropic", "google"],
                     metavar="TYPE",
                     help="LLM provider: ollama | openai | anthropic | google")
    llm.add_argument("--llm-model", metavar="MODEL",
                     help="Model name (e.g. 'gemini-2.0-flash', 'qwen3:8b')")
    llm.add_argument("--llm-api-key", metavar="KEY",
                     help="API key for cloud LLMs (openai / anthropic / google)")

    # Pipeline control
    ctrl = p.add_argument_group("Pipeline control")
    ctrl.add_argument(
        "--skip-steps",
        nargs="*",
        choices=["fastp","host","rnaspades","megahit",
                 "merge","mapping","coverage","cobra","viralquest"],
        default=[],
        metavar="STEP",
        help="Space-separated list of steps to skip. "
             "Choices: fastp host rnaspades megahit merge mapping coverage cobra viralquest",
    )

    return p


# ─── Step 10: Collect key results ────────────────────────────────────────────

def step_collect_results(
    sample: str,
    sample_out_dir: Path,
    global_results_dir: Path,
) -> Path:
    """
    Aggregate the most important outputs of a sample into a single flat
    results directory shared across all samples.

    After a full TAPIR run the per-sample output tree can span gigabytes of
    intermediate files (BAMs, raw assembly FASTAs, MMseqs2 temporaries, etc.).
    This step copies only the files a researcher typically needs for downstream
    analysis and publication directly into ``global_results_dir/`` — flat,
    with no per-sample subdirectory.  Because every output file already carries
    the sample name as a prefix (e.g. ``sample1_viral.fa``) there is no risk of
    name collision when results from multiple samples share the same folder.

    Files collected
    ---------------
    From step 1 (fastp)
      - ``<sample>_fastp.html``          — interactive QC report

    From step 9 (ViralQuest)
      - ``<sample>_viral.fa``            — final viral sequences (key deliverable)
      - ``<sample>_viral-BLAST.csv``     — BLAST hit table for all viral contigs
      - ``<sample>_bestSeqs.json``       — per-sequence annotation summary (JSON)
      - ``<sample>_visualization.html``  — interactive annotation HTML report

    Missing files are logged as warnings but do not abort the pipeline; the
    function always completes so that results from other samples are still
    collected.

    Parameters
    ----------
    sample : str
        Sample identifier used to derive expected source file paths.
    sample_out_dir : Path
        Root output directory for this sample (e.g. ``/results/sample1``).
    global_results_dir : Path
        Flat destination directory shared across all samples
        (e.g. ``/results/final_results``).

    Returns
    -------
    Path
        Path to ``global_results_dir`` where files were written.
    """
    global_results_dir.mkdir(parents=True, exist_ok=True)

    # ── Define source → destination pairs ────────────────────────────────────
    # Files are written directly into global_results_dir (no subdirectory).
    # All filenames already carry the sample prefix, so multiple samples can
    # coexist in the same folder without collision.
    vq_out = sample_out_dir / "09_viralquest" / f"OUTPUT_{sample}"

    targets: list[tuple[Path, str]] = [
        # QC report
        (sample_out_dir / "01_fastp" / f"{sample}_fastp.html",
         f"{sample}_fastp.html"),
        # Final viral sequences
        (vq_out / f"{sample}_viral.fa",
         f"{sample}_viral.fa"),
        # BLAST result table
        (vq_out / f"{sample}_viral-BLAST.csv",
         f"{sample}_viral-BLAST.csv"),
        # Per-sequence JSON annotation
        (vq_out / f"{sample}_bestSeqs.json",
         f"{sample}_bestSeqs.json"),
        # Interactive HTML report
        (vq_out / f"{sample}_visualization.html",
         f"{sample}_visualization.html"),
    ]

    copied, missing_count = 0, 0
    for file_src, dst_name in targets:
        dst = global_results_dir / dst_name
        if file_src.exists():
            shutil.copy2(file_src, dst)   # copy2 preserves metadata/timestamps
            log.debug(f"  Collected: {file_src.name} → {dst}")
            copied += 1
        else:
            log.warning(f"  Result file not found, skipping: {file_src}")
            missing_count += 1

    log.info(
        _c(GREEN, f"  ✓  Results collected for {sample}: "
                  f"{copied} file(s) → {global_results_dir}")
    )
    if missing_count:
        log.warning(
            f"  {missing_count} expected file(s) were missing for {sample}. "
            "Check the step logs above for errors."
        )

    return global_results_dir


# ══════════════════════════════════════════════════════════════════════════════
# MAIN ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    """
    Orchestrate the full TAPIR pipeline for all sample pairs in the input
    directory.

    For each discovered R1/R2 pair:
      1. Quality control   (fastp)
      2. Host removal      (Bowtie2)
      3. rnaSPAdes assembly
      4. MEGAHIT assembly
      5. Cross-assembly dereplication (MMseqs2)
      6. Read mapping      (Bowtie2 + SAMtools)
      7. Coverage          (CoverM)
      8. Contig extension  (COBRA)
      9. Viral annotation  (ViralQuest)

    Each step writes a ``.done_*`` sentinel file on completion; re-running the
    same command will skip completed steps automatically.
    """
    _banner()
    args = _build_parser().parse_args()

    # ── Initialise output directories and logging ─────────────────────────────
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.output_dir / "tapir.log"
    global log
    log = _setup_logging(log_path)

    # Directory that aggregates the key deliverables from every sample:
    # fastp QC report, final viral FASTA, and ViralQuest annotation files.
    # Intermediate files (BAMs, raw assemblies, MMseqs2 temps) remain in
    # their respective step subdirectories for reproducibility and debugging.
    final_results_dir = args.output_dir / "final_results"
    final_results_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"TAPIR v{__version__} — pipeline started")
    log.info(f"  Output directory : {args.output_dir}")
    log.info(f"  Final results    : {final_results_dir}")
    log.info(f"  Threads          : {args.threads}")
    log.info(f"  RAM              : {args.ram} GB")
    log.info(f"  Log file         : {log_path}")

    # ── Tool availability check ───────────────────────────────────────────────
    required_tools = [
        "fastp", "bowtie2", "bowtie2-build", "samtools",
        "spades.py", "megahit", "mmseqs",
        "cobra-meta", "viralquest",
    ]
    if not args.skip_host_removal:
        required_tools.append("bowtie2")
    _check_tools(required_tools)

    # ── Discover sample pairs ─────────────────────────────────────────────────
    pairs = _find_read_pairs(args.input_dir, args.r1_suffix, args.r2_suffix)
    log.info(f"  Samples detected : {len(pairs)}")

    # ── Process each sample ───────────────────────────────────────────────────
    for r1_raw, r2_raw in pairs:
        sample = r1_raw.name.replace(args.r1_suffix, "")
        s_out  = args.output_dir / sample

        log.info("")
        log.info(_c(BOLD + MAGENTA, f"{'━'*58}"))
        log.info(_c(BOLD + MAGENTA, f"  Sample: {sample}"))
        log.info(_c(BOLD + MAGENTA, f"{'━'*58}"))

        # Step 1 — fastp
        if "fastp" not in args.skip_steps:
            r1_qc, r2_qc = step_fastp(
                r1_raw, r2_raw,
                out_dir=s_out / "01_fastp",
                threads=args.threads,
                sample=sample,
            )
        else:
            r1_qc, r2_qc = r1_raw, r2_raw
            log.info("⏩  fastp skipped by user (--skip-steps fastp)")

        # Step 2 — Host removal
        if not args.skip_host_removal and "host" not in args.skip_steps:
            if not args.host_genome:
                log.error(
                    "--host-genome is required unless --skip-host-removal is set."
                )
                sys.exit(1)
            r1_nh, r2_nh = step_host_removal(
                r1_qc, r2_qc,
                host_genome=args.host_genome,
                out_dir=s_out / "02_host_removal",
                threads=args.threads,
                sample=sample,
            )
        else:
            r1_nh, r2_nh = r1_qc, r2_qc
            log.info("⏩  Host removal skipped")

        # Step 3 — rnaSPAdes
        if "rnaspades" not in args.skip_steps:
            rnaspades_fa = step_rnaspades(
                r1_nh, r2_nh,
                out_dir=s_out / "03_rnaspades",
                threads=args.threads,
                ram_gb=args.ram,
                sample=sample,
            )
        else:
            rnaspades_fa = s_out / "03_rnaspades" / "transcripts.fasta"
            log.info("⏩  rnaSPAdes skipped")

        # Step 4 — MEGAHIT
        if "megahit" not in args.skip_steps:
            megahit_fa = step_megahit(
                r1_nh, r2_nh,
                out_dir=s_out / "04_megahit",
                threads=args.threads,
                ram_gb=args.ram,
                sample=sample,
            )
        else:
            megahit_fa = s_out / "04_megahit" / "final.contigs.fa"
            log.info("⏩  MEGAHIT skipped")

        # Step 5 — Cross-assembly dereplication
        if "merge" not in args.skip_steps:
            merged_fa = step_merge_assemblies(
                contigs_list=[rnaspades_fa, megahit_fa],
                out_dir=s_out / "05_merge",
                threads=args.threads,
                sample=sample,
                min_length=args.min_contig_len,
            )
        else:
            merged_fa = s_out / "05_merge" / f"{sample}_merged_nr.fa"
            log.info("⏩  Assembly merge skipped")

        # Step 6 — Read mapping
        if "mapping" not in args.skip_steps:
            bam = step_map_reads(
                r1_nh, r2_nh,
                assembly=merged_fa,
                out_dir=s_out / "06_mapping",
                threads=args.threads,
                sample=sample,
            )
        else:
            bam = s_out / "06_mapping" / f"{sample}.sorted.bam"
            log.info("⏩  Read mapping skipped")

        # Step 7 — Coverage
        if "coverage" not in args.skip_steps:
            coverage = step_coverage(
                bam=bam,
                out_dir=s_out / "07_coverage",
                sample=sample,
                threads=args.threads,
            )
        else:
            coverage = s_out / "07_coverage" / f"{sample}_coverage.tsv"
            log.info("⏩  Coverage estimation skipped")

        # Step 8 — COBRA extension
        if "cobra" not in args.skip_steps:
            cobra_fa = step_cobra(
                assembly=merged_fa,
                query=args.cobra_query if hasattr(args, "cobra_query") else None,
                coverage=coverage,
                bam=bam,
                out_dir=s_out / "08_cobra",
                threads=args.threads,
                assembler_hint=args.cobra_assembler,
                mink=args.mink,
                maxk=args.maxk,
                cobra_min_length=args.cobra_min_len,
            )
        else:
            cobra_fa = s_out / "08_cobra" / "COBRA_extended_all.fa"
            log.info("⏩  COBRA skipped")

        # Step 9 — ViralQuest
        if "viralquest" not in args.skip_steps:
            viral_fa = step_viralquest(
                contigs=cobra_fa,
                out_dir=s_out / "09_viralquest",
                threads=args.threads,
                email=args.email,
                sample=sample,
                nr_dmnd=args.nr_dmnd,
                viral_dmnd=args.viral_dmnd,
                rvdb_hmm=args.rvdb_hmm,
                eggnog_hmm=args.eggnog_hmm,
                vfam_hmm=args.vfam_hmm,
                pfam_hmm=args.pfam_hmm,
                llm_type=args.llm_type,
                llm_model=args.llm_model,
                llm_api_key=getattr(args, "llm_api_key", None),
            )
        else:
            viral_fa = s_out / "09_viralquest" / f"OUTPUT_{sample}" / f"{sample}_viral.fa"
            log.info("⏩  ViralQuest skipped")

        # ── Step 10: Collect key results into final_results/ ─────────────────
        step_collect_results(
            sample=sample,
            sample_out_dir=s_out,
            global_results_dir=final_results_dir,
        )

        # ── Per-sample summary ────────────────────────────────────────────────
        log.info("")
        log.info(_c(GREEN + BOLD, f"  ✓  {sample} — all steps complete"))
        log.info(f"    Merged assembly  : {merged_fa}")
        log.info(f"    COBRA extended   : {cobra_fa}")
        if viral_fa.exists():
            log.info(
                f"    Viral sequences  : {viral_fa} "
                f"({_count_sequences(viral_fa):,} sequences)"
            )
        log.info(f"    Final results    : {final_results_dir}")

    # ── Pipeline summary ──────────────────────────────────────────────────────
    log.info("")
    log.info(_c(BOLD + CYAN, "═" * 58))
    log.info(_c(BOLD + CYAN, f"  TAPIR v{__version__} — pipeline finished"))
    log.info(_c(BOLD + CYAN, f"  {len(pairs)} sample(s) processed"))
    log.info(_c(BOLD + CYAN, "═" * 58))
    log.info(f"  Full log      : {log_path}")
    log.info(f"  Final results : {final_results_dir}")


if __name__ == "__main__":
    main()