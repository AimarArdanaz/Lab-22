#!/usr/bin/env python3
"""
count_genes.py — counts reads per gene with featureCounts over the 44 BAM files
mapped against F1, and builds the count matrix for DESeq2.

Why this and not the BigWig tracks
----------------------------------
The IGV profiles are there to LOOK at. To say "this gene is expressed more in
biofilm" you need one number per gene per sample, which is what featureCounts
produces, and then a statistical test (DESeq2) that takes the variability
between replicates into account.

Decisions this script makes, and why
------------------------------------
- Annotation: IGV/F1_anotacion.gff3 (the Bakta one, converted). CDS, tRNA,
  rRNA, ncRNA and tmRNA are counted; rRNA is counted ON PURPOSE, because
  knowing what fraction of each library the rRNA takes up is part of quality
  control (a sample that is 60 % rRNA holds far less useful information).
- Identifier: locus_tag (they are unique in this GFF3).
- Strandness (-s): NOT assumed. It is measured by measure_strandness(), which
  tries the three conventions on one sample and keeps whichever assigns the
  most reads.
- Paired-end: FRAGMENTS are counted (-p --countReadPairs). -B is not used
  (it demands both mates aligned properly) because in this data 44-50 % of the
  pairs came out discordant (bowtie2 ran without -X, maximum fragment 500 bp),
  and requiring it would throw away half the data.
- Single-end and paired-end are counted separately and then merged: in
  single-end each number is a read, in paired-end a fragment. It makes no
  difference to DESeq2 (it normalizes per column), but it is worth knowing.

Usage: python3 count_genes.py
Output in counts/:
   counts_single.tsv, counts_paired.tsv   raw featureCounts output
   counts_matrix.tsv                      genes x samples, ready for DESeq2
   assignment_summary.tsv                 % assigned and % rRNA per sample
"""
import os, subprocess, sys, csv, shutil

ROOT = os.path.dirname(os.path.abspath(__file__))
GFF = os.path.join(ROOT, "IGV", "F1_anotacion.gff3")   # file name as it is on disk
OUTDIR = os.path.join(ROOT, "counts")
TYPES = "CDS,tRNA,rRNA,ncRNA,tmRNA"
THREADS = "4"

if not shutil.which("featureCounts"):
    sys.exit("featureCounts is missing. Install it with:\n"
             "  conda activate tesis && conda install -c bioconda subread")

def is_paired(bam):
    """A BAM is paired if it has reads with flag 0x1 (read is paired)."""
    n = subprocess.run(["samtools", "view", "-c", "-f", "1", bam],
                       stdout=subprocess.PIPE, text=True, check=True).stdout
    return int(n.strip() or 0) > 0

def samples():
    for folder in ("biofilm", "no_biofilm"):
        d = os.path.join(ROOT, folder)
        for run in sorted(os.listdir(d)):
            bam = os.path.join(d, run, f"{run}_sorted.bam")
            if run.startswith("SRR") and os.path.exists(bam):
                yield run, bam

def featurecounts(bams, dest, strand, paired):
    cmd = ["featureCounts", "-a", GFF, "-F", "GTF", "-t", TYPES, "-g", "locus_tag",
           "-s", str(strand), "-T", THREADS, "-o", dest]
    if paired:
        cmd += ["-p", "--countReadPairs"]
    cmd += bams
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        raise RuntimeError("featureCounts failed")
    return r.stdout

def assigned(summary):
    """Reads featureCounts' .summary file and returns (assigned, total)."""
    data = {}
    with open(summary) as fh:
        header = fh.readline().split("\t")[1:]
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            data[fields[0]] = [int(x) for x in fields[1:]]
    asg = data["Assigned"]
    total = [sum(v[i] for v in data.values()) for i in range(len(asg))]
    return asg, total, [c.strip() for c in header]

def measure_strandness(bam, paired, label):
    """Tries -s 0/1/2 and keeps whichever assigns the most reads.
    0 = unstranded, 1 = forward (ligation), 2 = reverse (dUTP, the usual one)."""
    print(f"\nMeasuring the strand convention on {label} ({os.path.basename(bam)})")
    tmp = os.path.join(OUTDIR, "_strand_test.txt")
    best, results = None, {}
    for s in (0, 1, 2):
        featurecounts([bam], tmp, s, paired)
        asg, total, _ = assigned(tmp + ".summary")
        pct = asg[0] / total[0] * 100 if total[0] else 0
        results[s] = pct
        print(f"   -s {s}: {pct:5.1f} % assigned")
        if best is None or pct > results[best]:
            best = s
    for f in (tmp, tmp + ".summary"):
        if os.path.exists(f):
            os.unlink(f)
    # If 'unstranded' wins only narrowly, the library is stranded anyway: we
    # compare the two directions against each other to decide.
    if best == 0 and max(results[1], results[2]) > results[0] * 0.8:
        best = 1 if results[1] > results[2] else 2
        print(f"   -> the library IS stranded; using -s {best}")
    else:
        print(f"   -> using -s {best}")
    return best

os.makedirs(OUTDIR, exist_ok=True)
everything = list(samples())
single = [(r, b) for r, b in everything if not is_paired(b)]
paired = [(r, b) for r, b in everything if is_paired(b)]
print(f"{len(single)} single-end samples, {len(paired)} paired-end")

summary_rows = []
matrix = {}       # locus_tag -> {sample: count}
gene_info = {}    # locus_tag -> (chr, start, end, strand, length)

for group, group_samples in (("single", single), ("paired", paired)):
    if not group_samples:
        continue
    is_p = group == "paired"
    strand = measure_strandness(group_samples[0][1], is_p, group)
    dest = os.path.join(OUTDIR, f"counts_{group}.tsv")
    print(f"\nCounting {len(group_samples)} {group} samples with -s {strand}...")
    featurecounts([b for _, b in group_samples], dest, strand, is_p)
    asg, total, _ = assigned(dest + ".summary")

    with open(dest) as fh:
        fh.readline()                      # command line
        header = fh.readline().rstrip("\n").split("\t")
        runs = [os.path.basename(c).replace("_sorted.bam", "") for c in header[6:]]
        for line in fh:
            c = line.rstrip("\n").split("\t")
            gene = c[0]
            gene_info[gene] = (c[1], c[2], c[3], c[4], c[5])
            matrix.setdefault(gene, {})
            for run, val in zip(runs, c[6:]):
                matrix[gene][run] = int(val)
    for i, run in enumerate(runs):
        summary_rows.append([run, group, strand, total[i], asg[i],
                             f"{asg[i]/total[i]*100:.1f}" if total[i] else "0"])

# % rRNA per sample: how much of the library the ribosomal RNA takes up
rrna_genes = set()
with open(GFF) as fh:
    for l in fh:
        if l.startswith("#"):
            continue
        c = l.split("\t")
        if len(c) > 8 and c[2] == "rRNA":
            lt = [a for a in c[8].split(";") if a.startswith("locus_tag=")]
            if lt:
                rrna_genes.add(lt[0].split("=", 1)[1])

order = sorted(matrix)
sample_order = [f[0] for f in summary_rows]
with open(os.path.join(OUTDIR, "counts_matrix.tsv"), "w") as fh:
    fh.write("locus_tag\t" + "\t".join(sample_order) + "\n")
    for gene in order:
        fh.write(gene + "\t" + "\t".join(str(matrix[gene].get(m, 0))
                                         for m in sample_order) + "\n")

with open(os.path.join(OUTDIR, "assignment_summary.tsv"), "w") as fh:
    fh.write("run\tlayout\tstrand_s\ttotal\tassigned\tpct_assigned\tpct_rRNA\n")
    for row in summary_rows:
        run = row[0]
        gene_total = sum(matrix[g].get(run, 0) for g in matrix)
        rr = sum(matrix[g].get(run, 0) for g in rrna_genes if g in matrix)
        row.append(f"{rr/gene_total*100:.1f}" if gene_total else "0")
        fh.write("\t".join(str(x) for x in row) + "\n")

print(f"\nMatrix: {len(order)} genes x {len(sample_order)} samples")
print(f"Written to {OUTDIR}")
