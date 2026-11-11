# Exercise solutions

Try each exercise before reading its solution. Run commands from `~/learn` unless stated otherwise.

[Back to the course map](README.md)

---

## Module 1: The Terminal

### Exercise 1

```bash
cd ~/learn
mkdir scripts results
ls
```

### Exercise 2

```bash
grep -c ">" data/example.fasta        # 4 sequences
grep -c "bap" data/example.fasta      # 2 headers
grep -i "BAP" data/example.fasta
```

`-i` ignores upper and lower case, so `BAP` also finds `bap`.

### Exercise 3

```bash
cp data/example.fasta results/bap_sequences.fasta
ls data/strain0*.fasta                # check first
rm data/strain0*.fasta
ls -R
```

### Exercise 4

```bash
tail -n +2 data/strains.csv | cut -d',' -f4 | sort | uniq -c > results/esp_prevalence.txt
tail -n +2 data/strains.csv | cut -d',' -f2,4 | sort | uniq -c >> results/esp_prevalence.txt
cat results/esp_prevalence.txt
```

```
      5 0
      5 1
      2 clinical,0
      4 clinical,1
      3 commensal,0
      1 commensal,1
```

*esp* prevalence is 5/10 (50%): 4/6 in clinical and 1/4 in commensal strains.

```bash
grep ",1,1$" data/strains.csv | wc -l   # 3 strains carry both genes
```

### Exercise 5

Add at the end of the script:

```bash
echo "--- Clinical strains ---"
grep -c "clinical" "$table"
```

Without an argument, `$1` is empty, so `table` is empty. The script prints `=== Analysis of:  ===` and `tail` fails because there is no file called "". The script does not check its input; Lesson 7 fixes that.

### Exercise 6

```bash
rm -f results/contig_counts.tsv
for f in data/genome_*.fasta; do
  name=$(basename "$f" .fasta)
  n=$(grep -c ">" "$f")
  printf "%s\t%s\n" "$name" "$n" >> results/contig_counts.tsv
done
cat results/contig_counts.tsv
```

```
genome_1	2
genome_2	2
genome_3	2
```

An alternative is to redirect the whole loop once: `done > results/contig_counts.tsv`.

### Exercise 7

Add after the line that computes `total`:

```bash
if [ "$total" -lt 5 ]; then
  echo "Warning: fewer than 5 strains, prevalence estimates will be very imprecise" >&2
fi
```

`data/strains_hospital2.csv` has exactly 5 strains, so no warning. A 3-strain table triggers it.

### Exercise 8

1. Prevalence of *esp* (column 4) by origin:

```bash
awk -F',' 'NR > 1 {n[$2]++; pos[$2] += $4}
  END {for (o in n) printf "%s\t%d/%d\t%.1f%%\n", o, pos[o], n[o], 100 * pos[o] / n[o]}' data/strains.csv
```

```
clinical	4/6	66.7%
commensal	1/4	25.0%
```

The order of the lines may vary: `for (o in n)` does not guarantee an order.

2. Header prefix:

```bash
sed 's/^>/>cohortA_/' data/example.fasta > results/example_renamed.fasta
```

3. *esp* without *bap*:

```bash
awk -F',' 'NR > 1 && $4 == 1 && $3 == 0 {print $1}' data/strains.csv
```

```
C04
C09
```

---

## Module 2: Bioinformatics tools

### Exercise 9

```bash
conda activate bioinfo
seqkit version
blastn -version
abricate --version
datasets --version
conda env export --from-history > ~/learn/environment.yml
cat ~/learn/environment.yml
```

The file lists the environment name, the channels and the packages you asked for. Anyone can recreate the environment with `conda env create -f environment.yml`.

### Exercise 10

```bash
seqkit seq -w 0 data/multiline.fasta > results/linear.fasta
seqkit grep -n -r -p "bap" data/example.fasta
seqkit translate -T 11 data/example.fasta
```

`strain04_control` ends with `*`: its last codon, `TAA`, is a stop codon.

### Exercise 11

Follow Lesson 11 with your species and `--limit 5`, then:

```bash
seqkit stats -T genomes/*.fna
```

`num_seqs` is the number of contigs (in a complete genome: the chromosome plus any plasmids) and `sum_len` the total length in bp.

### Exercise 12

Make a FASTA file with only the Esp protein, then run the search loop from Lesson 12 with it:

```bash
seqkit grep -p esp references/amyloid_proteins.faa > references/esp.faa
```

```bash
rm -f results/esp_hits.tsv
for g in genomes/*.fna; do
  name=$(basename "$g" .fna)
  tblastn -query references/esp.faa -db "blastdb/${name}" \
    -outfmt "6 qseqid sseqid pident length qlen qcovs evalue bitscore" -evalue 1e-10 \
  | awk -v s="$name" 'BEGIN {OFS = "\t"} {print s, $0}' >> results/esp_hits.tsv
done
```

```bash
awk -F'\t' '$4 >= 90 && $7 >= 80 {print $1}' results/esp_hits.tsv | sort -u | wc -l
```

`sort -u` keeps each genome once, even if it has several hits.

### Exercise 13

The linking column is `strain`: the genome name in the presence matrix must match exactly the `strain` value in `metadata.csv`. Check for mismatches early, for example with `cut -d',' -f1` on both files.

---

## Module 3: R

### Exercise 14

```r
esp <- c(0, 1, 0, 1, 0, 1, 0, 0, 1, 1)
mean(esp) * 100                                     # 50

mean(strains$bap[strains$origin == "clinical"]) * 100   # 83.3

table(strains$bap, strains$esp)
```

The last table crosses *bap* (rows) with *esp* (columns): how many strains have neither, only one or both genes. Three strains carry both.

### Exercise 15

```r
strains <- strains |> mutate(both = as.integer(bap == 1 & esp == 1))

strains |>
  group_by(origin) |>
  summarise(n = n(), both = sum(both), prevalence = 100 * mean(both))
```

All 3 double-positive strains are clinical: 3/6 (50%) in clinical and 0/4 (0%) in commensal strains.

### Exercise 16

```r
binom.test(5, 6)$conf.int
```

The 95% CI goes from about 36% to 99.6%. With 6 strains, the true prevalence could be almost anything above one third.

```r
fisher.test(table(strains$origin, strains$bap))
```

p is about 0.19: not significant. 83% versus 25% looks like a big difference, but with 10 strains it could easily arise by chance. This is why sample size matters and why you report confidence intervals.

### Exercise 17

```r
prev_ci |>
  filter(gene == "esp") |>
  ggplot(aes(x = origin, y = prevalence)) +
  geom_col(width = 0.6, fill = "steelblue") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Origin", y = "esp prevalence (%)") +
  theme_classic()

ggsave("results/esp_prevalence.pdf", width = 5, height = 4)
```

---

## Module 4: Python

### Exercise 18

```python
def reverse_complement(seq):
    pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(pairs[base] for base in seq.upper()[::-1])

def count_codon(seq, codon):
    seq = seq.upper()
    codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
    return codons.count(codon.upper())

print(reverse_complement("ATGAAAGCG"))                  # CGCTTTCAT
print(count_codon("ATGAAAGCGAAA", "AAA"))              # 2
```

### Exercise 19

```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

proteins = []
for record in SeqIO.parse("data/example.fasta", "fasta"):
    record.seq = record.seq.translate(table=11)
    proteins.append(record)
SeqIO.write(proteins, "results/example_proteins.faa", "fasta")

high_gc = [r for r in SeqIO.parse("data/example.fasta", "fasta") if gc_fraction(r.seq) > 0.5]
SeqIO.write(high_gc, "results/high_gc.fasta", "fasta")
```

Two sequences pass the GC filter: `strain01_bap` (53.8%) and `strain03_bap` (56.4%).

### Exercise 20

With 90% identity and 90% coverage, `GCF_C` loses its *bap* call (its best hit has 85.3% identity):

| strain | bap (80/70) | esp (80/70) | bap (90/90) | esp (90/90) |
|---|---|---|---|---|
| GCF_A | 1 | 0 | 1 | 0 |
| GCF_B | 0 | 1 | 0 | 1 |
| GCF_C | 1 | 0 | 0 | 0 |
| GCF_D | 0 | 1 | 0 | 1 |
| GCF_E | 0 | 0 | 0 | 0 |

*bap* prevalence drops from 40% to 20%. Thresholds change results, so always report and justify them.

```python
print(matrix.sum(axis=1))
```

With 80/70, `GCF_A` to `GCF_D` carry one gene each and `GCF_E` none.

---

## Module 5: Good practice

### Exercise 21

```bash
cd ~/learn
git init
# create .gitignore as in Lesson 21
git add .gitignore scripts/
git commit -m "Add course scripts"
gh repo create learn-bioinformatics --private --source . --push
```

### Exercise 22

```bash
mkdir -p ~/thesis-amyloid/{data/raw,references,scripts,results,logs}
cd ~/thesis-amyloid
git init
```

Write `README.md` with the question, the data sources and the command to run (`bash run_all.sh`), then:

```bash
git add README.md
git commit -m "Create project structure"
```

Git does not track empty folders. Add an empty `.gitkeep` file to folders you want to keep: `touch scripts/.gitkeep`.
