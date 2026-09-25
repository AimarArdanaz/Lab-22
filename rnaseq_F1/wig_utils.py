#!/usr/bin/env python3
"""
wig_utils.py
====================================================================
Functions shared by pipeline_single_end.py and pipeline_paired_end.py
to split a sorted BAM into strands (+ / -) and convert the resulting
coverage into WIG format (variableStep).

This file must be in the SAME folder as those two scripts,
since both do `from wig_utils import ...`.

No new dependencies added: everything is done with `samtools`, which
you already had installed for the rest of the pipeline.

--------------------------------------------------------------------
WHAT IS "STRANDED COVERAGE" AND WHY IS IT DONE THIS WAY HERE?
--------------------------------------------------------------------
When you sequence RNA-seq with a "stranded" kit, each read retains
information about which DNA strand was originally transcribed.
Bowtie2 flags in the SAM/BAM FLAG whether a read mapped to the +
strand or the - strand (bit 0x10, value 16 of the flag).

BUT bit 16 of the flag is NOT directly "the transcribed strand" in
every protocol:

- dUTP kits / Illumina TruSeq Stranded mRNA (the MOST COMMON kits
  commercially today): the read you sequence is COMPLEMENTARY to
  the original template strand -> you have to INVERT the mapping
  strand to get the real transcribed strand.
  This is what we call 'fr-firststrand' here.

- Older "Ligation" kits (rare nowadays): the read DOES match the
  transcribed strand as mapped.
  This is what we call 'fr-secondstrand' here.

- If your library isn't "stranded" (standard genomic protocol, or
  non-directional RNA-seq), there's no strand to recover: use
  'unstranded' and a single WIG with total coverage is generated.

IF YOU DON'T KNOW WHICH OF THE TWO YOU USED (dUTP vs Ligation), the
most reliable way to check is to run the RSeQC package's
`infer_experiment.py` tool ONCE on one of your already-aligned BAMs,
using a BED annotation of your genome:

    pip install RSeQC
    infer_experiment.py -i SRRxxxx_sorted.bam -r annotation.bed

It returns the % of reads consistent with each convention: whichever
has the high percentage (usually >80-90%) is your real protocol.
For bacteria, you can generate the BED from a Prokka/NCBI GFF (there
are converters like BEDOPS's `gff2bed`).

Why does this matter for your thesis? Because in bacterial genomes
genes are very tightly packed on both strands, with frequent overlaps
and antisense arrangements -- exactly the case for pathogenicity
islands like the PAI of E. faecalis, where esp is located. Without
splitting by strand you can't tell whether the coverage signal
belongs to the gene of interest or to something overlapping on the
opposite strand.
--------------------------------------------------------------------
"""

import subprocess
from pathlib import Path


# ---------------------------------------------------------------------------
# Low-level samtools utilities
# ---------------------------------------------------------------------------

def _samtools_view_filter(input_bam, flag_f, flag_F, output_bam, threads):
    """
    Writes to output_bam the reads from input_bam that satisfy:
      - flag_f: bits that MUST be present (samtools -f), or None
      - flag_F: bits that must NOT be present (samtools -F), or None
    """
    cmd = ["samtools", "view", "-@", str(threads), "-b"]
    if flag_f is not None:
        cmd += ["-f", str(flag_f)]
    if flag_F is not None:
        cmd += ["-F", str(flag_F)]
    cmd += [str(input_bam)]
    with open(output_bam, "wb") as out:
        subprocess.run(cmd, stdout=out, check=True)


def _merge_bams(bam_list, output_bam, threads):
    """
    Merges 2+ already coordinate-sorted BAMs into one, preserving
    order (they all come from the same original sorted BAM, so the
    result stays sorted without re-sorting).
    """
    cmd = ["samtools", "merge", "-f", "-@", str(threads), str(output_bam)]
    cmd += [str(b) for b in bam_list]
    subprocess.run(cmd, check=True)


def compute_scale_factor(bam_path, threads=4):
    """
    Returns the scale factor to normalize to CPM
    (Counts Per Million mapped reads): 1e6 / total_mapped_reads.
    Computed on the whole BAM (not per strand), so that the + and -
    strands of the same sample end up on the same scale.
    """
    result = subprocess.run(
        ["samtools", "view", "-@", str(threads), "-c", "-F", "4", str(bam_path)],
        stdout=subprocess.PIPE, text=True, check=True,
    )
    total_mapped = int(result.stdout.strip() or 0)
    if total_mapped == 0:
        return 1.0
    return 1_000_000 / total_mapped


def bam_to_wig(bam_path, wig_path, threads=4, scale_factor=1.0, track_name=None):
    """
    Converts a (already filtered/sorted) BAM into a WIG (variableStep,
    span=1) with per-position coverage depth.

    Only positions with depth > 0 are written (default behavior of
    'samtools depth'), which is valid in WIG format: positions not
    listed are interpreted as having no data.

    If scale_factor != 1.0, values are normalized (e.g. to CPM)
    instead of writing raw depth.
    """
    name = track_name or Path(wig_path).stem
    proc = subprocess.Popen(
        ["samtools", "depth", "-@", str(threads), str(bam_path)],
        stdout=subprocess.PIPE,
        text=True,
    )

    current_chrom = None
    with open(wig_path, "w") as wig:
        wig.write(f'track type=wiggle_0 name="{name}"\n')
        for line in proc.stdout:
            chrom, pos, depth = line.rstrip("\n").split("\t")
            if chrom != current_chrom:
                wig.write(f"variableStep chrom={chrom} span=1\n")
                current_chrom = chrom
            if scale_factor != 1.0:
                wig.write(f"{pos} {float(depth) * scale_factor:.3f}\n")
            else:
                wig.write(f"{pos} {depth}\n")

    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(f"samtools depth failed for {bam_path}")


# ---------------------------------------------------------------------------
# Strand splitting: SINGLE-END
# ---------------------------------------------------------------------------

def split_strands_single_end(sorted_bam, prefix, threads, library_strand,
                              normalize_cpm=False):
    """
    Generates the coverage WIG(s) from an already sorted and indexed
    single-end BAM.

    For an individual read, the mapping strand is read directly from
    flag 16 (0x10 = aligned on the reverse strand):
      - flag 16 NOT present -> read mapped on strand +
      - flag 16 present     -> read mapped on strand -

    library_strand:
      'fr-firststrand'  -> dUTP / TruSeq Stranded kit (the most common):
                            the mapping strand is INVERTED relative
                            to the real transcribed strand.
      'fr-secondstrand' -> Ligation-type kit: the mapping strand
                            matches the transcribed strand.
      'unstranded'      -> no splitting; a single WIG with total coverage.

    Returns the list of Paths to the generated WIGs.
    """
    sorted_bam = Path(sorted_bam)
    outdir = sorted_bam.parent
    scale = compute_scale_factor(sorted_bam, threads) if normalize_cpm else 1.0

    if library_strand == "unstranded":
        wig = outdir / f"{prefix}.wig"
        bam_to_wig(sorted_bam, wig, threads, scale, track_name=prefix)
        return [wig]

    fwd_bam = outdir / f"{prefix}_fwdmapped.tmp.bam"  # reads mapped on strand +
    rev_bam = outdir / f"{prefix}_revmapped.tmp.bam"  # reads mapped on strand -

    # -F 20 = excludes unmapped(4) and reverse(16)  -> leaves only + mapped reads
    _samtools_view_filter(sorted_bam, None, 20, fwd_bam, threads)
    # -f 16 = requires reverse(16)                  -> leaves only - mapped reads
    _samtools_view_filter(sorted_bam, 16, None, rev_bam, threads)

    if library_strand == "fr-firststrand":
        # the mapping is inverted relative to the real transcribed strand
        plus_source, minus_source = rev_bam, fwd_bam
    else:  # fr-secondstrand
        plus_source, minus_source = fwd_bam, rev_bam

    plus_wig = outdir / f"{prefix}_plus.wig"
    minus_wig = outdir / f"{prefix}_minus.wig"
    bam_to_wig(plus_source, plus_wig, threads, scale, track_name=f"{prefix}_plus")
    bam_to_wig(minus_source, minus_wig, threads, scale, track_name=f"{prefix}_minus")

    fwd_bam.unlink(missing_ok=True)
    rev_bam.unlink(missing_ok=True)

    return [plus_wig, minus_wig]


# ---------------------------------------------------------------------------
# Strand splitting: PAIRED-END
# ---------------------------------------------------------------------------

def split_strands_paired_end(sorted_bam, prefix, threads, library_strand,
                              normalize_cpm=False):
    """
    Same as split_strands_single_end but for paired-end.

    In paired-end, the FRAGMENT's strand of origin isn't read from a
    single read: the flags of both mates need to be combined (this is
    what tools like RSeQC call the
    '1+-,1-+,2++,2--' / '1++,1--,2+-,2-+' convention):

      reads_A = (mate2 mapped on +) OR (mate1 mapped on -)
      reads_B = (mate2 mapped on -) OR (mate1 mapped on +)

    With 'fr-firststrand' (dUTP, the most common case), reads_A
    corresponds to the real transcribed + strand; with
    'fr-secondstrand' it's the other way round.

    SAM flags used:
      64  = mate1 (first in pair)
      128 = mate2 (second in pair)
      16  = aligned on the reverse strand
    """
    sorted_bam = Path(sorted_bam)
    outdir = sorted_bam.parent
    scale = compute_scale_factor(sorted_bam, threads) if normalize_cpm else 1.0

    if library_strand == "unstranded":
        wig = outdir / f"{prefix}.wig"
        bam_to_wig(sorted_bam, wig, threads, scale, track_name=prefix)
        return [wig]

    a1 = outdir / f"{prefix}_a1.tmp.bam"  # mate2 forward  (-f 128 -F 16)
    a2 = outdir / f"{prefix}_a2.tmp.bam"  # mate1 reverse  (-f 80,  80=64+16)
    b1 = outdir / f"{prefix}_b1.tmp.bam"  # mate2 reverse  (-f 144, 144=128+16)
    b2 = outdir / f"{prefix}_b2.tmp.bam"  # mate1 forward  (-f 64  -F 16)

    _samtools_view_filter(sorted_bam, 128, 16, a1, threads)
    _samtools_view_filter(sorted_bam, 80, None, a2, threads)
    _samtools_view_filter(sorted_bam, 144, None, b1, threads)
    _samtools_view_filter(sorted_bam, 64, 16, b2, threads)

    reads_A = outdir / f"{prefix}_readsA.tmp.bam"
    reads_B = outdir / f"{prefix}_readsB.tmp.bam"
    _merge_bams([a1, a2], reads_A, threads)
    _merge_bams([b1, b2], reads_B, threads)

    if library_strand == "fr-firststrand":
        plus_source, minus_source = reads_A, reads_B
    else:  # fr-secondstrand
        plus_source, minus_source = reads_B, reads_A

    plus_wig = outdir / f"{prefix}_plus.wig"
    minus_wig = outdir / f"{prefix}_minus.wig"
    bam_to_wig(plus_source, plus_wig, threads, scale, track_name=f"{prefix}_plus")
    bam_to_wig(minus_source, minus_wig, threads, scale, track_name=f"{prefix}_minus")

    for tmp in (a1, a2, b1, b2, reads_A, reads_B):
        tmp.unlink(missing_ok=True)

    return [plus_wig, minus_wig]
