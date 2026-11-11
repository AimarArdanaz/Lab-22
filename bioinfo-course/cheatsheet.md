# Cheat sheet

[Back to the course map](README.md)

## Terminal: navigation and files

| Command | What it does |
|---|---|
| `pwd` | Show the current folder |
| `ls`, `ls -l`, `ls -R` | List files; with details; recursively |
| `cd folder`, `cd ..`, `cd ~` | Enter a folder; go up; go home |
| `mkdir -p a/b/c` | Create folders, including parents |
| `touch file` | Create an empty file |
| `cp a b`, `cp -r dir1 dir2` | Copy a file; copy a folder |
| `mv a b`, `mv -i a b` | Move or rename; asking before overwriting |
| `rm -i file`, `rm -ri dir` | Delete with confirmation (no Bin!) |
| `open .` | Open the current folder in Finder |
| `Ctrl + C` | Stop the running command |

## Terminal: reading and processing text

| Command | What it does |
|---|---|
| `cat`, `less`, `head -n 5`, `tail -n 5` | Show a file; scroll; first or last lines |
| `wc -l` | Count lines |
| `grep "x" file` | Lines containing `x` |
| `grep -c`, `-v`, `-i`, `-A1`, `-q` | Count; invert; ignore case; plus next line; silent |
| `cut -d',' -f2,3` | Keep columns 2 and 3 of a CSV |
| `sort`, `sort -k2,2n`, `uniq -c` | Sort; by column 2 numerically; count repeats |
| `tr ',' '\t'` | Replace characters (CSV to TSV) |
| `sed 's/a/b/g'` | Replace text |
| `sed -i '' 's/a/b/' file` | Edit in place (macOS) |
| `awk -F',' 'NR > 1 {print $1}'` | Print column 1, skipping the header |

## Redirection and pipes

| Symbol | What it does |
|---|---|
| `>` | Save output to a file (overwrite) |
| `>>` | Append output to a file |
| `<` | Use a file as input |
| `<< 'EOF'` | Use the following lines as input, until `EOF` |
| `\|` | Send output to the next command |
| `2>&1` | Include error messages in the output |
| `\| tee file` | Show on screen and save to a file |

## Scripts

```bash
#!/bin/bash
set -euo pipefail

input=$1                                  # first argument

if [ ! -f "$input" ]; then                # check it exists
  echo "Error: $input not found" >&2
  exit 1
fi

for f in data/*.fasta; do                 # loop over files
  name=$(basename "$f" .fasta)
  echo "$name: $(grep -c '>' "$f")"
done

while read -r id; do                      # loop over lines of a file
  echo "$id"
done < ids.txt
```

## conda

| Command | What it does |
|---|---|
| `conda create -n env tool1 tool2` | Create an environment |
| `conda activate env` / `conda deactivate` | Enter / leave it |
| `conda install tool` | Install in the active environment |
| `conda env list` | List environments |
| `conda env export --from-history > environment.yml` | Save the recipe |
| `conda env create -f environment.yml` | Recreate from the recipe |

## seqkit

| Command | What it does |
|---|---|
| `seqkit stats -a -T *.fasta` | Summary statistics |
| `seqkit fx2tab -n -l -g in.fa` | Name, length, GC per sequence |
| `seqkit seq -w 0 in.fa` | One line per sequence |
| `seqkit seq -m 500 in.fa` | Keep sequences of at least 500 bp |
| `seqkit grep -n -r -p "bap" in.fa` | Sequences whose header matches |
| `seqkit translate -T 11 in.fa` | Translate with the bacterial code |

## BLAST and abricate

```bash
makeblastdb -in genome.fna -dbtype nucl -out blastdb/genome

tblastn -query proteins.faa -db blastdb/genome \
  -outfmt "6 qseqid sseqid pident length qlen qcovs evalue bitscore" \
  -evalue 1e-10

abricate --datadir abricate_db --setupdb
abricate --datadir abricate_db --db amyloid --minid 80 --mincov 60 genomes/*.fna > hits.tab
abricate --summary hits.tab > summary.tab
```

## R

```r
library(tidyverse)
df <- read_csv("data/strains.csv")

df |> filter(origin == "clinical")
df |> group_by(origin) |> summarise(n = n(), prev = 100 * mean(bap))
df |> pivot_longer(c(bap, esp), names_to = "gene", values_to = "present")
left_join(meta, matrix, by = "strain")

binom.test(5, 6)$conf.int                 # exact 95% CI
fisher.test(table(df$origin, df$bap))     # compare groups
p.adjust(p_values, method = "BH")         # multiple testing

ggplot(df, aes(x, y, fill = group)) + geom_col() + theme_classic()
ggsave("figure.pdf", width = 6, height = 4)
```

## Python

```python
from Bio import SeqIO
import pandas as pd

for r in SeqIO.parse("in.fasta", "fasta"):
    print(r.id, len(r.seq), r.seq.translate(table=11))

df = pd.read_csv("file.tsv", sep="\t", names=cols)
df[(df["pident"] >= 80) & (df["qcovs"] >= 70)]
df.pivot_table(index="strain", columns="gene", values="present", aggfunc="max", fill_value=0)
df.groupby("origin")[["bap", "esp"]].mean() * 100
df.to_csv("out.csv", index=False)
```

## Git

| Command | What it does |
|---|---|
| `git init` | Start tracking a folder |
| `git status` | See what changed |
| `git add file` | Stage a change |
| `git commit -m "message"` | Save a snapshot |
| `git log --oneline` | See the history |
| `git push` / `git pull` | Send to / bring from GitHub |
| `gh repo create name --private --source . --push` | Create a GitHub repo and push |
