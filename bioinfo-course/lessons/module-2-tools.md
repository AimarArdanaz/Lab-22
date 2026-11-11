# Module 2: Bioinformatics tools

Now you install real bioinformatics software, download real genomes and search them for the genes of interest. By the end you will have a presence/absence matrix: one row per genome, one column per gene.

[Back to the course map](../README.md)

- [Lesson 9: conda environments](#lesson-9-conda-environments)
- [Lesson 10: seqkit](#lesson-10-seqkit)
- [Lesson 11: Downloading genomes from NCBI](#lesson-11-downloading-genomes-from-ncbi)
- [Lesson 12: Local BLAST](#lesson-12-local-blast)
- [Lesson 13: Abricate with a custom database](#lesson-13-abricate-with-a-custom-database)

---

## Lesson 9: conda environments

**Goal:** install tools in an isolated environment you can recreate on any computer.

### Why environments?

Each project may need different versions of the same tool. An **environment** is a separate folder with its own programs. The name of the active environment appears in brackets at the start of the prompt, for example `(base)` or `(bioinfo)`.

### Create and use an environment

```bash
conda create -n bioinfo seqkit blast abricate ncbi-datasets-cli
```

Type `y` when asked. If you already have an environment called `bioinfo`, install into it instead:

```bash
conda activate bioinfo
```

```bash
conda install seqkit blast abricate ncbi-datasets-cli
```

Check that the tools work:

```bash
seqkit version
```

```bash
blastn -version
```

```bash
abricate --version
```

```bash
datasets --version
```

### Everyday commands

| Command | What it does |
|---|---|
| `conda activate bioinfo` | Enter the environment |
| `conda deactivate` | Leave it |
| `conda env list` | List your environments |
| `conda list` | List the packages in the active environment |
| `conda install name` | Install a package |
| `conda env export --from-history > environment.yml` | Save the recipe of the environment |
| `conda env create -f environment.yml` | Recreate the environment from the recipe |

The `environment.yml` file is what makes your analysis reproducible: commit it with your scripts (Lesson 21).

### Apple Silicon note

Most Bioconda tools now exist for Apple Silicon (`osx-arm64`). If conda says `PackagesNotFoundError` for a tool, create an Intel environment, which runs through Rosetta:

```bash
softwareupdate --install-rosetta
```

```bash
CONDA_SUBDIR=osx-64 conda create -n bioinfo-x86 abricate blast
```

```bash
conda activate bioinfo-x86
```

```bash
conda config --env --set subdir osx-64
```

The last command makes future installs in that environment also use Intel packages.

### Exercise 9

Activate your environment, check the four tools above, and export the environment to `~/learn/environment.yml`. Open the file with `cat`: what does it contain?

---

## Lesson 10: seqkit

**Goal:** inspect, filter and transform FASTA files correctly, including multi-line sequences.

Real FASTA files split each sequence across many lines. Create one:

```bash
cat > data/multiline.fasta << 'EOF'
>seqA strain01 putative surface protein
ATGAAAGCGTTACCGATTGCG
CTGGCAGTGACCGGTATG
>seqB strain02 hypothetical protein
ATGAGTAACGTTCTGAAAGCC
GGTCTGATCGCAGCGTTA
EOF
```

`grep -A1` from Lesson 4 would now return half sequences. `seqkit` handles this properly.

### Summary statistics

```bash
seqkit stats data/*.fasta
```

```bash
seqkit stats -a -T data/*.fasta > results/fasta_stats.tsv
```

`-a` adds extra columns (N50, GC content...), `-T` gives a tab-separated table.

### Per-sequence table

```bash
seqkit fx2tab -n -l -g data/multiline.fasta
```

`-n` name only (no sequence), `-l` length, `-g` GC content.

### Common operations

| Task | Command |
|---|---|
| One line per sequence | `seqkit seq -w 0 in.fasta` |
| Keep sequences of at least 500 bp | `seqkit seq -m 500 in.fasta` |
| Sequences whose header matches a pattern | `seqkit grep -n -r -p "bap" in.fasta` |
| Sequences whose ID is in a list | `seqkit grep -f ids.txt in.fasta` |
| Find a motif | `seqkit locate -p ATGAAAGCG in.fasta` |
| Extract a region | `seqkit subseq -r 1:12 in.fasta` |
| Translate to protein | `seqkit translate -T 11 in.fasta` |
| Reverse complement | `seqkit seq -r -p in.fasta` |
| Remove duplicate sequences | `seqkit rmdup -s in.fasta` |
| Sort by length, longest first | `seqkit sort -l -r in.fasta` |
| Add a prefix to headers | `seqkit replace -p '^' -r 'strain01_' in.fasta` |

In `seqkit grep`, `-n` matches the full header, `-r` treats the pattern as a regular expression and `-p` gives the pattern. `-T 11` uses the **bacterial genetic code** (table 11).

### Exercise 10

1. Convert `data/multiline.fasta` to one line per sequence and save it as `results/linear.fasta`.
2. Extract the *bap* sequences from `data/example.fasta` with `seqkit grep`.
3. Translate `data/example.fasta` with the bacterial code. Which sequence ends with a stop codon (`*`)?

---

## Lesson 11: Downloading genomes from NCBI

**Goal:** download a set of genome assemblies and organise them with clean names.

The `datasets` tool from NCBI downloads genomes by species, taxon or accession. As an example we use *Enterococcus faecalis*, the species in which *esp* (enterococcal surface protein) was described. For *bap* (biofilm-associated protein) the classic species is *Staphylococcus aureus*. Use the species of your own study.

### 1. Explore what is available

```bash
datasets summary genome taxon "Enterococcus faecalis" --assembly-source RefSeq --assembly-level complete --limit 10 --as-json-lines | dataformat tsv genome --fields accession,organism-name,assminfo-level
```

`dataformat` turns NCBI's JSON output into a table. To see every available field:

```bash
dataformat tsv genome --help
```

### 2. Save a list of accessions

```bash
mkdir -p ~/learn/project && cd ~/learn/project
```

```bash
datasets summary genome taxon "Enterococcus faecalis" --assembly-source RefSeq --assembly-level complete --limit 10 --as-json-lines | dataformat tsv genome --fields accession --elide-header > accessions.txt
```

```bash
cat accessions.txt
```

`--elide-header` removes the column name, leaving a clean list.

### 3. Download them

```bash
datasets download genome accession --inputfile accessions.txt --include genome --filename genomes.zip
```

```bash
unzip -o genomes.zip -d ncbi_download
```

The genomes end up in `ncbi_download/ncbi_dataset/data/<accession>/<accession>_<name>_genomic.fna`. Copy them into one folder, named by accession:

```bash
mkdir -p genomes
for f in ncbi_download/ncbi_dataset/data/GC*/*_genomic.fna; do
  acc=$(basename "$(dirname "$f")")
  cp "$f" "genomes/${acc}.fna"
done
```

```bash
ls genomes
```

```bash
seqkit stats genomes/*.fna
```

### Your own isolates

If you sequenced your own clinical strains, put their assemblies in the same `genomes/` folder with the `.fna` extension and a short name without spaces, such as `HUN_0012.fna`.

### Metadata

Isolation source, host or country live in each genome's **BioSample** record and are often incomplete. For clinical versus commensal comparisons you will usually build `metadata.csv` yourself (from BioSample, the original papers or your hospital records) with at least these columns:

```
strain,origin
GCF_000000000.1,clinical
```

**Never commit identifiable patient data to GitHub.** Use anonymous codes.

### Exercise 11

Download 5 complete genomes of your species of interest, rename them by accession and use `seqkit stats` to answer: how many contigs (or chromosomes plus plasmids) does each one have, and what is its total length?

---

## Lesson 12: Local BLAST

**Goal:** search every genome for the reference proteins and keep only convincing hits.

### Reference sequences

Create `references/amyloid_proteins.faa` with the protein sequences of the genes you study, in FASTA format, downloaded from UniProt or NCBI Protein. Use short headers with the gene name:

```
>bap
MKSKRLLAAALTAG...
>esp
MNKKWFLSAAL...
```

(The sequences above are truncated placeholders.) Write down the accession of every reference in your methods. Bap and Esp belong to a family of large surface proteins with repeat regions, so homologues in other species can differ in length and sequence: choose references that fit your species.

### Why `tblastn`?

| Program | Query | Database |
|---|---|---|
| `blastn` | DNA | DNA |
| `blastp` | Protein | Protein |
| `tblastn` | **Protein** | **DNA translated in 6 frames** |

`tblastn` finds genes in raw genome assemblies even when the DNA has diverged, because proteins change more slowly than DNA.

### 1. Build one database per genome

```bash
mkdir -p blastdb
for g in genomes/*.fna; do
  name=$(basename "$g" .fna)
  makeblastdb -in "$g" -dbtype nucl -out "blastdb/${name}"
done
```

### 2. Search

```bash
mkdir -p results
rm -f results/blast_hits.tsv

for g in genomes/*.fna; do
  name=$(basename "$g" .fna)
  tblastn -query references/amyloid_proteins.faa \
    -db "blastdb/${name}" \
    -outfmt "6 qseqid sseqid pident length qlen qcovs evalue bitscore" \
    -evalue 1e-10 \
    -max_target_seqs 5 \
  | awk -v s="$name" 'BEGIN {OFS = "\t"} {print s, $0}' >> results/blast_hits.tsv
done
```

- `\` at the end of a line continues the command on the next line.
- `rm -f` empties the results first, because the loop appends with `>>`.
- `awk -v s="$name"` passes the genome name into `awk`, which adds it as the first column.

### 3. Understand the columns

| # | Column | Meaning |
|---|---|---|
| 1 | strain | Genome name (added by `awk`) |
| 2 | qseqid | Reference protein (`bap`, `esp`) |
| 3 | sseqid | Contig where the hit is |
| 4 | pident | Percent identity |
| 5 | length | Alignment length |
| 6 | qlen | Length of the reference protein |
| 7 | qcovs | Percent of the reference covered by all hits on that contig |
| 8 | evalue | Expected number of chance hits |
| 9 | bitscore | Score; higher is better |

### 4. Filter and keep the best hit

Keep hits with at least 80% identity and 70% coverage:

```bash
awk -F'\t' '$4 >= 80 && $7 >= 70' results/blast_hits.tsv > results/blast_filtered.tsv
```

Keep the best hit per genome and gene (sort by genome, gene and bitscore, then keep the first line of each pair):

```bash
sort -t$'\t' -k1,1 -k2,2 -k9,9gr results/blast_filtered.tsv | awk -F'\t' '!seen[$1 FS $2]++' > results/blast_best.tsv
```

**About thresholds:** there is no universal cut-off. Genes with long repeat regions, like *bap* and *esp*, may be split across contigs in short-read assemblies, which lowers coverage. Look at the distribution of `pident` and `qcovs` before choosing, and justify the values in your methods.

### Exercise 12

Run `tblastn` with only the Esp reference on your downloaded genomes. How many genomes have a hit with at least 90% identity and 80% coverage?

---

## Lesson 13: Abricate with a custom database

**Goal:** screen all genomes in one command and get a presence/absence matrix directly.

`abricate` is a wrapper around BLAST designed for exactly this task. It comes with resistance and virulence databases, and you can add your own. It searches **nucleotide** sequences, so it complements the protein-level search of Lesson 12.

### 1. Prepare the database

Abricate needs nucleotide sequences with headers in this format:

```
>DATABASE~~~GENE~~~ACCESSION~~~DESCRIPTION
```

For example (sequences truncated):

```
>amyloid~~~bap~~~ACCESSION_1~~~Biofilm-associated protein
ATGAAAAAG...
>amyloid~~~esp~~~ACCESSION_2~~~Enterococcal surface protein
ATGAATAAA...
```

Save this file as `references/amyloid_genes.fna`. Keep the database inside your project folder so it travels with your scripts:

```bash
mkdir -p abricate_db/amyloid
```

```bash
cp references/amyloid_genes.fna abricate_db/amyloid/sequences
```

```bash
abricate --datadir abricate_db --setupdb
```

```bash
abricate --datadir abricate_db --list
```

The file must be called exactly `sequences`.

### 2. Screen the genomes

```bash
abricate --datadir abricate_db --db amyloid --minid 80 --mincov 60 genomes/*.fna > results/abricate_amyloid.tab
```

### 3. Summary matrix

```bash
abricate --summary results/abricate_amyloid.tab > results/amyloid_summary.tab
```

```bash
cat results/amyloid_summary.tab
```

Each row is a genome; each gene column holds the coverage (%) when found and `.` when absent.

### 4. Convert to 0/1

```bash
awk -F'\t' 'BEGIN {OFS = ","}
  NR == 1 {$1 = "strain"; print; next}
  {sub(/.*\//, "", $1); sub(/\.fna$/, "", $1)
   for (i = 3; i <= NF; i++) $i = ($i == "." ? 0 : 1)
   print}' results/amyloid_summary.tab > results/presence_matrix.csv
```

- `sub(/.*\//, "", $1)` removes the folder from the file name, and the second `sub` removes `.fna`.
- `($i == "." ? 0 : 1)` means "if the cell is `.` write 0, otherwise 1".
- Setting `OFS = ","` writes the result as CSV.

The second column, `NUM_FOUND`, is the number of genes found in each genome.

### Abricate or tblastn?

Use both and compare. Abricate (DNA) is fast and standard; `tblastn` (protein) catches divergent alleles that abricate misses. Report which one you used for the final calls, and the thresholds.

### Exercise 13

Join `results/presence_matrix.csv` with your `metadata.csv` by strain. You can do it in R (Lesson 15) or Python (Lesson 20), so write down for now which column links the two files.

Next: [Module 3, R](module-3-r.md)
