#!/usr/bin/env python3
"""
Pipeline SRA -> FASTQ -> Alignment (PAIRED-END) -> sorted BAM -> WIG
====================================================================
Usage:
    python pipeline_paired_end.py srr_list.txt

Requirements (must be in PATH):
    - sra-tools (prefetch, fasterq-dump) >= 3.x
    - bowtie2
    - samtools

IMPORTANT: this script imports wig_utils.py, which must be in the
same folder (it contains the strand-splitting and WIG-generation
functions; it adds no new dependencies, everything uses samtools).

The srr_list.txt file must have one SRR accession per line, e.g.:
    SRR38278097
    SRR38278098
    SRR38278099
"""

import subprocess
import sys
import shutil
from pathlib import Path
from datetime import datetime

from wig_utils import split_strands_paired_end

# ---------------------- CONFIGURATION (edit here) ------------------------
REF_INDEX = "../F1_index"   # path to the bowtie2 index (without the .bt2 extension)
THREADS = 4                 # number of threads for fasterq-dump / bowtie2 / samtools
OUTDIR = Path(".")          # output folder for FASTQ/SAM/BAM/WIG
KEEP_FASTQ = False          # True = keep the .fastq files after aligning
KEEP_SRA = False            # True = keep the folder with the .sra after extracting

# --- Stranded coverage configuration ----------------------------------------
# 'fr-firststrand'  -> dUTP kit / Illumina TruSeq Stranded (the MOST COMMON
#                       commercially available today): the strand inferred
#                       from the mates' flags is INVERTED relative to the
#                       actually transcribed strand.
# 'fr-secondstrand' -> Ligation-type kit (less common today): matches
#                       the transcribed strand directly.
# 'unstranded'      -> don't split by strand; a single WIG with total coverage.
#
# If you don't know which one you used, check your sequencing provider's
# kit datasheet, or check empirically with RSeQC (see wig_utils.py,
# infer_experiment.py function, on one of your already-generated BAMs).
LIBRARY_STRAND = "fr-firststrand"

# True = normalizes coverage to CPM (counts per million mapped reads)
# instead of raw depth. Useful if you're going to compare the WIG of one
# sample against another with different sequencing depth.
NORMALIZE_CPM = False
# -----------------------------------------------------------------------------


def log(msg, logfile):
    """Prints to screen and saves to the log file."""
    print(msg)
    with open(logfile, "a") as f:
        f.write(msg + "\n")


def run(cmd, logfile, capture_to=None):
    """
    Runs an external command. If capture_to is given, redirects stdout+stderr
    to that file. Returns the exit code (returncode).
    """
    log(f"$ {' '.join(cmd)}", logfile)
    if capture_to:
        with open(capture_to, "w") as out:
            result = subprocess.run(cmd, stdout=out, stderr=out)
    else:
        with open(logfile, "a") as out:
            result = subprocess.run(cmd, stdout=out, stderr=out)
    return result.returncode


def main():
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <srr_list_file.txt>")
        sys.exit(1)

    lista_path = Path(sys.argv[1])
    if not lista_path.is_file():
        print(f"ERROR: file '{lista_path}' not found")
        sys.exit(1)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    logfile = OUTDIR / f"pipeline_log_{timestamp}.txt"

    log(f"=== Paired-end pipeline started: {datetime.now()} ===", logfile)
    log(f"    Strand mode (LIBRARY_STRAND): {LIBRARY_STRAND}", logfile)
    log(f"    CPM normalization: {NORMALIZE_CPM}", logfile)

    # Read the list, stripping \r (in case the .txt comes from Windows) and blank lines
    with open(lista_path, "r") as f:
        sra_list = [line.strip().replace("\r", "") for line in f if line.strip()]

    total = len(sra_list)
    ok = 0
    fail = 0

    for i, srr in enumerate(sra_list, start=1):
        log("", logfile)
        log("=" * 60, logfile)
        log(f">>> [{i}/{total}] Processing {srr}", logfile)
        log("=" * 60, logfile)

        # ---------- 1. Download (prefetch) ----------
        run(["prefetch", srr, "-O", str(OUTDIR)], logfile)
        sra_file = OUTDIR / srr / f"{srr}.sra"
        if not sra_file.is_file():
            log(f"!!! ERROR: prefetch failed for {srr}, skipping.", logfile)
            fail += 1
            continue

        # ---------- 2. SRA -> FASTQ (split into _1 / _2) ----------
        run(
            [
                "fasterq-dump", str(sra_file),
                "--split-files", "-e", str(THREADS), "-O", str(OUTDIR),
            ],
            logfile,
        )
        fastq_1 = OUTDIR / f"{srr}_1.fastq"
        fastq_2 = OUTDIR / f"{srr}_2.fastq"

        if not fastq_1.is_file() or not fastq_2.is_file():
            log(
                f"!!! ERROR: fasterq-dump did not generate both mates for {srr} "
                f"(is it really paired-end?), skipping.",
                logfile,
            )
            fail += 1
            continue

        # ---------- 3. Alignment with bowtie2 (paired-end) ----------
        sam = OUTDIR / f"{srr}.sam"
        summary = OUTDIR / f"{srr}_bowtie2_summary.txt"

        run(
            [
                "bowtie2", "-p", str(THREADS), "-x", REF_INDEX,
                "-1", str(fastq_1),
                "-2", str(fastq_2),
                "-S", str(sam),
            ],
            logfile,
            capture_to=summary,
        )

        if not sam.is_file() or sam.stat().st_size == 0:
            log(f"!!! ERROR: bowtie2 did not generate a SAM for {srr}, skipping.", logfile)
            fail += 1
            continue

        log(f"Bowtie2 summary ({srr}):", logfile)
        with open(summary) as f:
            log(f.read(), logfile)

        # ---------- 4. SAM -> sorted, indexed BAM ----------
        bam = OUTDIR / f"{srr}_sorted.bam"

        view_proc = subprocess.Popen(
            ["samtools", "view", "-@", str(THREADS), "-b", "-F", "4", str(sam)],
            stdout=subprocess.PIPE,
        )
        subprocess.run(
            ["samtools", "sort", "-@", str(THREADS), "-o", str(bam)],
            stdin=view_proc.stdout,
        )
        view_proc.stdout.close()
        view_proc.wait()

        run(["samtools", "index", str(bam)], logfile)

        # ---------- 5. BAM -> WIG (coverage, split by strand) ----------
        try:
            wig_files = split_strands_paired_end(
                bam, srr, THREADS, LIBRARY_STRAND, normalize_cpm=NORMALIZE_CPM
            )
            for w in wig_files:
                log(f"    WIG generated -> {w}", logfile)
        except Exception as e:
            log(f"!!! ERROR generating WIG for {srr}: {e}", logfile)
            # Not marked as a total failure: the BAM was generated correctly.

        # ---------- 6. Cleanup ----------
        sam.unlink(missing_ok=True)

        if not KEEP_FASTQ:
            fastq_1.unlink(missing_ok=True)
            fastq_2.unlink(missing_ok=True)

        if not KEEP_SRA:
            shutil.rmtree(OUTDIR / srr, ignore_errors=True)

        log(f">>> {srr} completed -> {bam}", logfile)
        ok += 1

    log("", logfile)
    log(f"=== Pipeline finished: {datetime.now()} ===", logfile)
    log(f"Completed successfully: {ok} / {total}", logfile)
    log(f"Failed: {fail} / {total}", logfile)


if __name__ == "__main__":
    main()
