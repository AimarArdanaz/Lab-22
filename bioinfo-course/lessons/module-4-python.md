# Module 4: Python for sequences and tables

Python is the most common general-purpose language in bioinformatics. It is ideal for working with sequences (Biopython) and for turning tool outputs into tables (pandas).

[Back to the course map](../README.md)

- [Lesson 18: Python basics](#lesson-18-python-basics)
- [Lesson 19: Biopython](#lesson-19-biopython)
- [Lesson 20: pandas](#lesson-20-pandas)

---

## Setup

Install Python and the libraries in your conda environment:

```bash
conda activate bioinfo
```

```bash
conda install python ipython biopython pandas
```

Three ways to run Python:

| Way | Command | Use it for |
|---|---|---|
| Interactive | `ipython` (quit with `exit`) | Trying things out |
| Script | `python scripts/my_script.py` | Real analyses |
| VS Code | Open a `.py` file and press the play button | Writing and running scripts |

Run everything from `~/learn` so relative paths work.

---

## Lesson 18: Python basics

**Goal:** learn the building blocks: variables, lists, dictionaries, loops, conditions and functions.

### Variables and types

```python
gene = "bap"            # text (str)
n_strains = 10          # integer (int)
prevalence = 60.0       # decimal (float)
is_present = True       # boolean (bool)

print(f"{gene}: {prevalence}% of {n_strains} strains")
```

An `f"..."` string inserts variables inside braces.

### Lists

```python
genes = ["bap", "esp", "agg"]
print(genes[0])         # first element: Python counts from 0
print(len(genes))
genes.append("ebpA")
```

### Dictionaries

A dictionary links keys to values:

```python
prevalence = {"bap": 60.0, "esp": 50.0}
print(prevalence["bap"])
prevalence["agg"] = 20.0

for gene, value in prevalence.items():
    print(gene, value)
```

### Loops and conditions

Indentation (4 spaces) defines the blocks. There are no `do`/`done` as in bash.

```python
for gene in genes:
    if gene == "bap":
        print(gene, "-> Staphylococcus")
    elif gene == "esp":
        print(gene, "-> Enterococcus")
    else:
        print(gene, "-> other")
```

### Functions

```python
def gc_content(seq):
    """Return the GC content of a DNA sequence, in percent."""
    seq = seq.upper()
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)

print(gc_content("ATGAAAGCGTTACCG"))
```

### Reading a FASTA file by hand

Understanding this helps you understand what Biopython does for you:

```python
sequences = {}
name = None

with open("data/example.fasta") as fh:
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            sequences[name] = ""
        else:
            sequences[name] += line

for name, seq in sequences.items():
    print(name, len(seq), round(gc_content(seq), 1))
```

`with open(...)` closes the file automatically. `strip()` removes the line break.

### A script with arguments

```python
# scripts/gc.py
import sys

def gc_content(seq):
    seq = seq.upper()
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)

seq = sys.argv[1]
print(f"GC = {gc_content(seq):.1f}%")
```

```bash
python scripts/gc.py ATGAAAGCGTTACCG
```

`sys.argv[1]` is the first argument, like `$1` in bash.

### Exercise 18

1. Write a function `reverse_complement(seq)` using a dictionary `{"A": "T", "T": "A", "G": "C", "C": "G"}`. Hint: `seq[::-1]` reverses a text.
2. Write a function that counts how many times a codon appears in a sequence, reading in frame (positions 0, 3, 6...).

---

## Lesson 19: Biopython

**Goal:** read, analyse, translate and write sequences with a library built for it.

### Read a FASTA file

```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

for record in SeqIO.parse("data/example.fasta", "fasta"):
    print(record.id, len(record.seq), round(100 * gc_fraction(record.seq), 1))
```

Each `record` has `.id`, `.description` and `.seq`. `SeqIO` handles multi-line sequences and many formats (`"genbank"`, `"fastq"`...).

### Work with sequences

```python
record = next(SeqIO.parse("data/example.fasta", "fasta"))
seq = record.seq

print(seq.reverse_complement())
print(seq.translate(table=11))
print(seq.translate(table=11, to_stop=True))
```

`table=11` is the bacterial genetic code. `to_stop=True` stops at the first stop codon.

### Filter and write

```python
bap_records = [r for r in SeqIO.parse("data/example.fasta", "fasta") if "bap" in r.id]
SeqIO.write(bap_records, "results/bap_biopython.fasta", "fasta")
```

The line `[r for r in ... if ...]` is a **list comprehension**: a compact loop that builds a list.

### A summary table of a genome

```python
import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

fasta = sys.argv[1]
print("contig\tlength\tgc")
for record in SeqIO.parse(fasta, "fasta"):
    print(f"{record.id}\t{len(record.seq)}\t{100 * gc_fraction(record.seq):.1f}")
```

Save it as `scripts/contig_table.py` and run it:

```bash
python scripts/contig_table.py data/example.fasta > results/contigs.tsv
```

### Exercise 19

1. Translate every sequence in `data/example.fasta` with table 11 and save the proteins to `results/example_proteins.faa`.
2. Write only the sequences with a GC content above 50% to `results/high_gc.fasta`.

---

## Lesson 20: pandas

**Goal:** turn BLAST hits into a presence/absence matrix and compute prevalences.

`data/example_blast_hits.tsv` contains made-up hits in the format of Lesson 12, for five genomes (`GCF_A` to `GCF_E`). `GCF_E` has no hits at all.

### Load and filter

```python
import pandas as pd

cols = ["strain", "qseqid", "sseqid", "pident", "length", "qlen", "qcovs", "evalue", "bitscore"]
hits = pd.read_csv("data/example_blast_hits.tsv", sep="\t", names=cols)

print(hits.head())
print(hits.shape)

good = hits[(hits["pident"] >= 80) & (hits["qcovs"] >= 70)]
```

`names=cols` gives column names because the file has no header. Conditions are combined with `&` (and) or `|` (or), each in parentheses.

### Best hit per strain and gene

```python
best = (good.sort_values("bitscore", ascending=False)
            .drop_duplicates(["strain", "qseqid"]))
```

### Presence/absence matrix

```python
matrix = (best.assign(present=1)
              .pivot_table(index="strain", columns="qseqid", values="present",
                           aggfunc="max", fill_value=0))
print(matrix)
```

**Trap:** genomes with no hits (`GCF_E`) and genes with no hits are missing from the matrix, and they would silently disappear from the prevalence. Add them back explicitly:

```python
all_strains = ["GCF_A", "GCF_B", "GCF_C", "GCF_D", "GCF_E"]
all_genes = ["bap", "esp"]

matrix = matrix.reindex(index=all_strains, columns=all_genes, fill_value=0)
print(matrix)
```

With real data, read the list of genomes from the folder:

```python
from pathlib import Path
all_strains = sorted(p.stem for p in Path("genomes").glob("*.fna"))
```

### Join metadata and compute prevalence

```python
meta = pd.read_csv("data/example_metadata.csv")
data = meta.merge(matrix, left_on="strain", right_index=True, how="left")

print(data)
print(data[["bap", "esp"]].mean() * 100)
print(data.groupby("origin")[["bap", "esp"]].mean() * 100)

data.to_csv("results/presence_with_metadata.csv", index=False)
```

The resulting CSV goes straight into R for statistics and plots (Lessons 16 and 17).

### Exercise 20

1. Build the matrix with stricter thresholds (90% identity, 90% coverage). What changes?
2. Count, for each strain, how many of the genes it carries (`matrix.sum(axis=1)`).

Next: [Module 5, good practice](module-5-good-practice.md)
