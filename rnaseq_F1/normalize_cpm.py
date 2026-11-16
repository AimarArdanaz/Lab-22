#!/usr/bin/env python3
"""
normalize_cpm.py — converts the raw-coverage WIG files into BigWig tracks
normalized to CPM (counts per million), so that track heights can be compared
BETWEEN samples.

The problem it solves
---------------------
The WIG/BigWig files produced by pipeline_single_end.py and
pipeline_paired_end.py are raw coverage: how many reads cover each base. That
depends on how deeply each sample was sequenced, and here the depth ranges from
1.3 to 49 million mapped reads (a 38x difference). Without correcting for it, a
track twice as tall as another one may be so only because it was sequenced
twice as deeply.

What it does
------------
It multiplies every coverage value by 1e6 / (mapped reads of THAT sample). The
result is "coverage per million mapped reads". The factor is computed on the
whole BAM (both strands together), exactly like compute_scale_factor() in
wig_utils.py, so that the + and - strands of the same sample end up on the same
scale.

What it does NOT fix (important for the thesis)
-----------------------------------------------
- Mappability bias between batches: SRR2915 and SRR3248 align at 95-98 %
  against F1 and the rest at 82-91 %. CPM does not correct that.
- Library composition: if a sample is dominated by rRNA, the CPM of every other
  gene goes down even though its real expression did not change. For serious
  quantification, count per gene and use DESeq2 normalization.
- CPM is for LOOKING at comparable profiles, not for calling differential
  expression.

It touches nothing that already exists: it reads the .wig files and writes new
.bw files into bw_cpm/.

Usage: python3 normalize_cpm.py            (all samples)
       python3 normalize_cpm.py SRR13200548 SRR32481956   (only those)
"""
import os, subprocess, sys, shutil, tempfile

ROOT = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(ROOT, "bw_cpm")
CHROM = os.path.join(ROOT, "F1.chrom.sizes")

if not shutil.which("wigToBigWig"):
    sys.exit("wigToBigWig is missing. Install it with:\n"
             "  conda activate tesis && conda install -c bioconda ucsc-wigtobigwig")

def samples():
    for folder in ("biofilm", "no_biofilm"):
        d = os.path.join(ROOT, folder)
        if not os.path.isdir(d):
            continue
        for run in sorted(os.listdir(d)):
            if run.startswith("SRR") and os.path.isdir(os.path.join(d, run)):
                yield run, os.path.join(d, run)

def cpm_factor(bam):
    """1e6 / mapped reads. idxstats reads the .bai index, so it is instant."""
    out = subprocess.run(["samtools", "idxstats", bam],
                         stdout=subprocess.PIPE, text=True, check=True).stdout
    mapped = sum(int(l.split("\t")[2]) for l in out.strip().split("\n") if l)
    if mapped == 0:
        raise RuntimeError(f"0 mapped reads in {bam}")
    return 1_000_000 / mapped, mapped

def scale(wig_in, factor, dest_bw):
    """Multiplies every value in the WIG by the factor and converts to BigWig."""
    with tempfile.NamedTemporaryFile("w", suffix=".wig", delete=False) as tmp:
        tmp_path = tmp.name
        with open(wig_in) as fh:
            for line in fh:
                if line[0].isdigit():
                    pos, val = line.split()
                    tmp.write(f"{pos} {float(val) * factor:.4f}\n")
                else:
                    # 'track ...' and 'variableStep ...' headers are copied as they are
                    tmp.write(line)
    try:
        subprocess.run(["wigToBigWig", tmp_path, CHROM, dest_bw], check=True)
    finally:
        os.unlink(tmp_path)

requested = set(sys.argv[1:])
os.makedirs(OUTDIR, exist_ok=True)
done = 0
for run, folder in samples():
    if requested and run not in requested:
        continue
    bam = os.path.join(folder, f"{run}_sorted.bam")
    f, mapped = cpm_factor(bam)
    for strand in ("plus", "minus"):
        wig = os.path.join(folder, f"{run}_{strand}.wig")
        bw = os.path.join(OUTDIR, f"{run}_{strand}_cpm.bw")
        scale(wig, f, bw)
    done += 1
    print(f"{run}: {mapped:,} mapped, factor {f:.6f} -> bw_cpm/{run}_*_cpm.bw",
          flush=True)
print(f"\n{done} samples normalized in {OUTDIR}")
