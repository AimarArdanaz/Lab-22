# Module 5: Good practice

Code that only works once, on one computer, is not enough for a thesis. This module covers version control with Git and a project layout that lets anyone (including you, in a year) rerun the whole analysis.

[Back to the course map](../README.md)

- [Lesson 21: Git and GitHub](#lesson-21-git-and-github)
- [Lesson 22: A reproducible project](#lesson-22-a-reproducible-project)

---

## Lesson 21: Git and GitHub

**Goal:** keep a history of every change to your scripts and back them up on GitHub.

### Concepts

| Term | Meaning |
|---|---|
| Repository (repo) | A folder whose history Git tracks |
| Commit | A saved snapshot, with a message describing the change |
| Staging | Choosing which changes go into the next commit (`git add`) |
| Remote | A copy of the repository on a server, such as GitHub |
| Push / pull | Send commits to the remote / bring them from it |

### One-time configuration

```bash
git config --global user.name "Your Name"
```

```bash
git config --global user.email "you@example.com"
```

```bash
git config --global init.defaultBranch main
```

### Daily workflow

```bash
cd ~/learn
```

```bash
git init
```

```bash
git status
```

```bash
git add scripts/prevalence.sh
```

```bash
git commit -m "Add prevalence script"
```

```bash
git log --oneline
```

Write commit messages that say **what changed and why**: "Filter BLAST hits at 80% identity", not "changes".

### What not to commit: `.gitignore`

Git is for code and small text files, not for genomes or large results. Create a `.gitignore` file listing what Git must ignore:

```bash
cat > .gitignore << 'EOF'
# Large or regenerable files
genomes/
ncbi_download/
blastdb/
*.zip
*.fna

# Results can be regenerated from the scripts
results/

# System files
.DS_Store
.Rhistory
.RData
.Rproj.user/
EOF
```

**Never commit patient data, identifiable metadata or passwords.** Once something is pushed to GitHub, removing it from the history is hard.

### Publish on GitHub

1. Create an account at https://github.com.
2. Install the GitHub command line tool and log in:

```bash
conda install gh
```

```bash
gh auth login
```

Choose GitHub.com, HTTPS, and log in with the browser.

3. Create the repository on GitHub and push in one step:

```bash
gh repo create my-thesis-analysis --private --source . --push
```

Use `--private` for thesis work until you decide to publish.

After that, each time you commit:

```bash
git push
```

### Exercise 21

Turn `~/learn` into a Git repository with a `.gitignore`, make your first commit with the scripts folder, and push it to a private GitHub repository.

---

## Lesson 22: A reproducible project

**Goal:** organise the thesis analysis so it can be rerun from raw data to figures with one command.

### Folder structure

```
thesis-amyloid/
├── README.md              # What the project is and how to run it
├── environment.yml        # Exact software (conda)
├── run_all.sh             # Runs every step in order
├── data/
│   ├── raw/               # Original data. Never edit by hand
│   └── metadata.csv
├── references/            # Reference sequences, with their accessions
├── scripts/
│   ├── 01_download_genomes.sh
│   ├── 02_search_genes.sh
│   ├── 03_build_matrix.py
│   └── 04_prevalence_stats.R
├── results/               # Everything here can be regenerated
└── logs/                  # Output of each run
```

### Principles

1. **Raw data is read-only.** Scripts read from `data/raw/` and write to `results/`. Protect raw files with `chmod -R a-w data/raw`.
2. **Number the scripts** in the order they run.
3. **One script, one task.** Short scripts are easier to check and rerun.
4. **Record the software:** `conda env export --from-history > environment.yml`, plus the version of each tool in your methods (`blastn -version`, `abricate --version`).
5. **Record the parameters:** thresholds (identity, coverage, e-value) at the top of the scripts as variables, not buried in commands.
6. **Keep logs.** `tee` shows the output on screen and saves it at the same time; `2>&1` includes error messages:

```bash
bash scripts/02_search_genes.sh 2>&1 | tee logs/02_search_genes.log
```

7. **Write a README** that explains how to run the project from scratch.

### Example `run_all.sh`

```bash
#!/bin/bash
set -euo pipefail

mkdir -p results logs

bash scripts/01_download_genomes.sh data/raw/accessions.txt 2>&1 | tee logs/01.log
bash scripts/02_search_genes.sh 2>&1 | tee logs/02.log
python scripts/03_build_matrix.py 2>&1 | tee logs/03.log
Rscript scripts/04_prevalence_stats.R 2>&1 | tee logs/04.log

echo "Pipeline finished"
```

`Rscript` runs an R script from the Terminal.

### Next step: workflow managers

When the project grows (hundreds of genomes, many steps), tools like **Snakemake** or **Nextflow** rerun only the steps whose inputs changed and run genomes in parallel. Learn them once the shell version of your pipeline works.

### Exercise 22

Create the folder structure above for your thesis with `mkdir -p`, add a `README.md` with a one-paragraph description and the steps to run, and commit it.

Next: [Final project](final-project.md)
