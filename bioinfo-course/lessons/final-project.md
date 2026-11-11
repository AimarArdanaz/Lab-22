# Final project: the full thesis pipeline

Put everything together: from a list of genomes to a figure and a table of prevalences with statistics, rerunnable with one command.

[Back to the course map](../README.md)

## The question

What is the prevalence of genes encoding amyloid-forming surface proteins (for example `bap` and `esp`) in genomes of clinical strains, and does it differ from commensal strains?

## Steps

| Step | Script | Input | Output | Lessons |
|---|---|---|---|---|
| 1 | `01_download_genomes.sh` | `data/raw/accessions.txt` | `genomes/*.fna` | 6, 11 |
| 2 | `02_search_genes.sh` | `genomes/`, `references/` | BLAST and abricate results | 7, 12, 13 |
| 3 | `03_build_matrix.py` | BLAST hits, `metadata.csv` | `results/presence_with_metadata.csv` | 20 |
| 4 | `04_prevalence_stats.R` | Presence table | Prevalence table, tests, figures | 15, 16, 17 |

Set up the project as in [Lesson 22](module-5-good-practice.md#lesson-22-a-reproducible-project).

## Step 1: download genomes

```bash
#!/bin/bash
# scripts/01_download_genomes.sh
# Usage: bash scripts/01_download_genomes.sh data/raw/accessions.txt
set -euo pipefail

if [ $# -ne 1 ] || [ ! -f "$1" ]; then
  echo "Usage: $0 <accessions.txt>" >&2
  exit 1
fi

accessions=$1
mkdir -p genomes

datasets download genome accession --inputfile "$accessions" --include genome --filename genomes.zip
unzip -o -q genomes.zip -d ncbi_download

for f in ncbi_download/ncbi_dataset/data/GC*/*_genomic.fna; do
  acc=$(basename "$(dirname "$f")")
  cp "$f" "genomes/${acc}.fna"
done

echo "Genomes ready: $(ls genomes/*.fna | wc -l | tr -d ' ')"
```

Your own assemblies go into `genomes/` as well.

## Step 2: search the genes

```bash
#!/bin/bash
# scripts/02_search_genes.sh
set -euo pipefail

# Parameters: write them in your methods
EVALUE=1e-10
MIN_ID=80
MIN_COV=60

mkdir -p blastdb results
rm -f results/blast_hits.tsv

# tblastn (protein level)
for g in genomes/*.fna; do
  name=$(basename "$g" .fna)
  makeblastdb -in "$g" -dbtype nucl -out "blastdb/${name}" > /dev/null
  tblastn -query references/amyloid_proteins.faa \
    -db "blastdb/${name}" \
    -outfmt "6 qseqid sseqid pident length qlen qcovs evalue bitscore" \
    -evalue "$EVALUE" -max_target_seqs 5 \
  | awk -v s="$name" 'BEGIN {OFS = "\t"} {print s, $0}' >> results/blast_hits.tsv
done

# abricate (nucleotide level)
abricate --datadir abricate_db --db amyloid --minid "$MIN_ID" --mincov "$MIN_COV" \
  genomes/*.fna > results/abricate_amyloid.tab
abricate --summary results/abricate_amyloid.tab > results/amyloid_summary.tab

echo "Search finished"
```

`> /dev/null` hides the messages of `makeblastdb`.

## Step 3: build the matrix

```python
# scripts/03_build_matrix.py
from pathlib import Path
import pandas as pd

MIN_PIDENT = 80
MIN_QCOVS = 70

cols = ["strain", "qseqid", "sseqid", "pident", "length", "qlen", "qcovs", "evalue", "bitscore"]
hits = pd.read_csv("results/blast_hits.tsv", sep="\t", names=cols)

good = hits[(hits["pident"] >= MIN_PIDENT) & (hits["qcovs"] >= MIN_QCOVS)]
best = good.sort_values("bitscore", ascending=False).drop_duplicates(["strain", "qseqid"])

all_strains = sorted(p.stem for p in Path("genomes").glob("*.fna"))
all_genes = [line[1:].split()[0] for line in open("references/amyloid_proteins.faa") if line.startswith(">")]

matrix = (best.assign(present=1)
              .pivot_table(index="strain", columns="qseqid", values="present",
                           aggfunc="max", fill_value=0)
              .reindex(index=all_strains, columns=all_genes, fill_value=0))

meta = pd.read_csv("data/metadata.csv")
data = meta.merge(matrix, left_on="strain", right_index=True, how="left")

missing = data[all_genes].isna().any(axis=1).sum()
if missing:
    print(f"Warning: {missing} strains in metadata have no genome")

data.to_csv("results/presence_with_metadata.csv", index=False)
print(data.groupby("origin")[all_genes].mean() * 100)
```

## Step 4: statistics and figures

```r
# scripts/04_prevalence_stats.R
library(tidyverse)

data <- read_csv("results/presence_with_metadata.csv")
genes <- setdiff(names(data), c("strain", "origin"))

long <- data |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "present") |>
  filter(!is.na(present))

prev <- long |>
  group_by(origin, gene) |>
  summarise(n = n(), positives = sum(present), prevalence = 100 * positives / n,
            .groups = "drop") |>
  rowwise() |>
  mutate(ci_low  = 100 * binom.test(positives, n)$conf.int[1],
         ci_high = 100 * binom.test(positives, n)$conf.int[2]) |>
  ungroup()

tests <- long |>
  group_by(gene) |>
  filter(sum(present) > 0, sum(present) < n()) |>
  summarise(p_value = fisher.test(table(origin, present))$p.value) |>
  mutate(p_adj = p.adjust(p_value, method = "BH"))

write_csv(prev, "results/prevalence_table.csv")
write_csv(tests, "results/fisher_tests.csv")

p <- ggplot(prev, aes(x = gene, y = prevalence, fill = origin)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Gene", y = "Prevalence (%)", fill = "Origin") +
  theme_classic()

ggsave("results/prevalence.pdf", p, width = 6, height = 4)
```

If your metadata has columns other than `strain` and `origin`, add them to `setdiff()`.

## Deliverables checklist

- [ ] `README.md` explaining the question and how to run the pipeline.
- [ ] `environment.yml` with the software.
- [ ] Reference sequences with their accessions listed in the README.
- [ ] `run_all.sh` runs from raw data to figures without manual steps.
- [ ] Prevalence table with n, positives, %, and 95% CI per group.
- [ ] Fisher tests with adjusted p-values.
- [ ] Figure saved as PDF.
- [ ] Comparison of tblastn and abricate calls, with the discrepancies explained.
- [ ] Everything committed to a private GitHub repository, with no patient data.

## Going further

- Add sequence type (MLST) to the metadata and check whether gene presence depends on lineage.
- Extract the gene sequences found and build a phylogeny of the alleles.
- Check whether the genes sit on plasmids or on the chromosome.
- Move the pipeline to Snakemake.
