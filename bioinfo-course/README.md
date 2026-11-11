# Bioinformatics from Scratch: from the Terminal to Gene Prevalence

A hands-on course for life scientists who have never programmed. It runs on **macOS** and goes step by step from opening the Terminal to a full, reproducible analysis: **how prevalent are genes encoding amyloid-forming surface proteins (such as `bap` and `esp`) in genomes of clinical bacterial strains?**

Each lesson has a goal, commands to copy and run, and an exercise. Solutions are in [`solutions.md`](solutions.md).

## Who this is for

- Biologists, biotechnologists and microbiologists starting with bioinformatics.
- No previous programming experience needed.
- You need a Mac (Apple Silicon or Intel) and an internet connection.

## Course map

| # | Module | Lesson | What you learn |
|---|---|---|---|
| 0 | Setup | [Prepare your Mac](lessons/module-0-setup.md) | Terminal, conda, a code editor |
| 1 | Terminal | [Navigating](lessons/module-1-terminal.md#lesson-1-navigating-the-file-system) | `pwd`, `ls`, `cd`, `mkdir` |
| 2 | Terminal | [Creating and reading files](lessons/module-1-terminal.md#lesson-2-creating-and-reading-files) | `cat`, `head`, `wc`, `grep` on a FASTA file |
| 3 | Terminal | [Copy, move, delete](lessons/module-1-terminal.md#lesson-3-copying-moving-and-deleting) | `cp`, `mv`, `rm`, wildcards |
| 4 | Terminal | [Redirection and pipes](lessons/module-1-terminal.md#lesson-4-redirection-and-pipes) | `>`, `>>`, `\|`, `cut`, `sort`, `uniq` |
| 5 | Terminal | [Your first script](lessons/module-1-terminal.md#lesson-5-your-first-script) | Variables, arguments, `chmod` |
| 6 | Terminal | [Loops](lessons/module-1-terminal.md#lesson-6-for-loops) | `for`, `while read` |
| 7 | Terminal | [Robust scripts](lessons/module-1-terminal.md#lesson-7-conditionals-and-robust-scripts) | `if`, exit codes, `set -euo pipefail` |
| 8 | Terminal | [awk and sed](lessons/module-1-terminal.md#lesson-8-awk-and-sed) | Processing tables and FASTA files |
| 9 | Tools | [conda environments](lessons/module-2-tools.md#lesson-9-conda-environments) | Installing bioinformatics software |
| 10 | Tools | [seqkit](lessons/module-2-tools.md#lesson-10-seqkit) | Inspecting and filtering sequences |
| 11 | Tools | [Downloading genomes](lessons/module-2-tools.md#lesson-11-downloading-genomes-from-ncbi) | NCBI `datasets` |
| 12 | Tools | [Local BLAST](lessons/module-2-tools.md#lesson-12-local-blast) | `makeblastdb`, `tblastn`, filtering hits |
| 13 | Tools | [Abricate](lessons/module-2-tools.md#lesson-13-abricate-with-a-custom-database) | Custom gene database, presence/absence matrix |
| 14 | R | [R basics](lessons/module-3-r.md#lesson-14-r-basics) | Vectors, data frames, reading CSV files |
| 15 | R | [tidyverse](lessons/module-3-r.md#lesson-15-prevalence-with-the-tidyverse) | Prevalence by group |
| 16 | R | [Statistics](lessons/module-3-r.md#lesson-16-statistics-for-prevalence-data) | Confidence intervals, Fisher's test, multiple testing |
| 17 | R | [ggplot2](lessons/module-3-r.md#lesson-17-plots-with-ggplot2) | Bar charts and heatmaps for publication |
| 18 | Python | [Python basics](lessons/module-4-python.md#lesson-18-python-basics) | Variables, lists, dictionaries, functions |
| 19 | Python | [Biopython](lessons/module-4-python.md#lesson-19-biopython) | Parsing and translating sequences |
| 20 | Python | [pandas](lessons/module-4-python.md#lesson-20-pandas) | From BLAST hits to a presence/absence matrix |
| 21 | Good practice | [Git and GitHub](lessons/module-5-good-practice.md#lesson-21-git-and-github) | Version control |
| 22 | Good practice | [Reproducible projects](lessons/module-5-good-practice.md#lesson-22-a-reproducible-project) | Folder structure, environments, logs |
| FP | All | [Final project](lessons/final-project.md) | The full thesis pipeline |

Also available: a [command cheat sheet](cheatsheet.md) and the [exercise solutions](solutions.md).

## How to use this course

1. Do the lessons in order. Each one builds on the previous one.
2. Type the commands yourself rather than only reading them. Copy and paste is fine at first.
3. Work inside one folder, `~/learn`, so every path in the lessons works as written.
4. When something fails, read the error message carefully. It usually says exactly what went wrong.

## Repository structure

```
.
├── README.md            # This page
├── cheatsheet.md        # Quick reference of commands
├── solutions.md         # Solutions to the exercises
├── data/                # Small example files used in the lessons
│   ├── example.fasta
│   ├── strains.csv
│   ├── example_blast_hits.tsv
│   └── example_metadata.csv
└── lessons/
    ├── module-0-setup.md
    ├── module-1-terminal.md
    ├── module-2-tools.md
    ├── module-3-r.md
    ├── module-4-python.md
    ├── module-5-good-practice.md
    └── final-project.md
```

## About the example data

All files in `data/` are **made up** for practice. The sequences are short invented fragments and the strain table is fictional. Real genomes are downloaded in [Lesson 11](lessons/module-2-tools.md#lesson-11-downloading-genomes-from-ncbi).
