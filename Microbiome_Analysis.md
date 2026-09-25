# 16S rRNA Analysis Pipeline — Synthetic Community

Complete, reproducible code for the analysis of the synthetic microbial community (9 species, based on Lekbua et al. 2024, *Cell Rep Methods* 4:100832), comparing three experimental conditions: **pHD17** (plasmid-carrying strain, functional *bap* gene), **Dbap** (plasmid-carrying strain, mutated *bap*) and **PF** (plasmid-free strain). 4 replicates per condition, 12 samples total, 16S rRNA sequencing.

This document collects, in execution order, all the code used from building the reference database through to the final differential abundance analysis. It is published alongside the manuscript as reproducible supplementary material.

## Table of contents

1. [Environment and dependencies](#1-environment-and-dependencies)
2. [Building the 16S reference database](#2-building-the-16s-reference-database)
3. [Primer trimming (cutadapt)](#3-primer-trimming-cutadapt)
4. [DADA2: from raw reads to ASVs](#4-dada2-from-raw-reads-to-asvs)
5. [Taxonomic assignment (SILVA + BLAST)](#5-taxonomic-assignment-silva--blast)
6. [Phylogenetic tree (optional)](#6-phylogenetic-tree-optional)
7. [Building the phyloseq object](#7-building-the-phyloseq-object)
8. [Quality and prevalence filtering](#8-quality-and-prevalence-filtering)
9. [Alpha diversity](#9-alpha-diversity)
10. [Beta diversity (Aitchison/CLR)](#10-beta-diversity-aitchisonclr)
11. [Rarefaction curves](#11-rarefaction-curves)
12. [Composition: bar plots, bubble plot, heatmap, Venn, UpSet](#12-composition-bar-plots-bubble-plot-heatmap-venn-upset)
13. [Differential abundance: LEfSe, ANOVA, ANCOM, ANCOM-BC2](#13-differential-abundance-lefse-anova-ancom-ancom-bc2)
14. [Known package issues and applied fixes](#14-known-package-issues-and-applied-fixes)
15. [Reproducibility](#15-reproducibility)

---

## 1. Environment and dependencies

### Command-line tools (conda/mamba)

```bash
mamba install -c conda-forge -c bioconda ncbi-datasets-cli entrez-direct barrnap seqkit unzip
mamba install -c bioconda cutadapt
mamba activate bioinfo
```

### R packages used

```r
library(dada2)
library(phyloseq)
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(ape)
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(plyr)
library(tidyverse)
library(data.table)
library(vegan)
library(zCompositions)
library(compositions)
library(microbiome)
library(microbiomeMarker)
library(ANCOMBC)
library(ampvis2)
library(UpSetR)
library(ggpubr)
library(see)
library(gghalves)
library(readxl)
library(knitr)
library(DT)
```

---

## 2. Building the 16S reference database

Before species-level taxonomy can be assigned, a curated reference database with the 16S gene of the 9 specific SkinCom strains is needed (generic SILVA is not enough — it only resolves well down to genus). `skincom_16S.sh` downloads each strain's genome from NCBI, extracts the 16S with `barrnap`, and filters by length (1300–1700 bp) to discard partial fragments.

```bash
#!/usr/bin/env bash
# =====================================================================
# skincom_16S.sh
# Retrieves the 16S rRNA of the 9 SkinCom strains (Lekbua et al. 2024,
# Cell Rep Methods 4:100832) from their genomes on NCBI.
#
# Dependencies (e.g. with conda/mamba):
#   mamba install -c conda-forge -c bioconda ncbi-datasets-cli \
#         entrez-direct barrnap seqkit unzip
#
# Usage:
#   bash skincom_16S.sh                # everything with default values
#   THREADS=8 OUT=my_folder bash skincom_16S.sh
#
# If a strain isn't found or the chosen genome isn't the right one,
# edit OUT/cepas.tsv (4th column = manual accession, e.g. GCF_000008345.1)
# and re-run the script: it reuses whatever was already downloaded.
# =====================================================================
set -uo pipefail

THREADS=${THREADS:-4}
OUT=${OUT:-skincom_16S}
MINLEN=1300     # the full 16S is ~1500 bp; discards partial hits
MAXLEN=1700

for x in datasets esearch efetch esummary xtract barrnap seqkit unzip awk; do
  command -v "$x" >/dev/null 2>&1 || { echo "ERROR: missing '$x' in PATH" >&2; exit 1; }
done

mkdir -p "$OUT"/{genomes,per_strain,fallback_nuccore,sin_match}

# ---------------------------------------------------------------------
# Strain table: label | species | strain terms (pipe-separated) | manual accession ("-" = automatic)
# ---------------------------------------------------------------------
if [[ ! -s "$OUT/cepas.tsv" ]]; then
  printf '%s\t%s\t%s\t%s\n' \
    "C_afermentans"   "Corynebacterium afermentans" "ATCC 51403|CIP 103499"      "-" \
    "C_acnes"         "Cutibacterium acnes"         "KPA171202|DSM 16379"         "-" \
    "M_luteus"        "Micrococcus luteus"          "ATCC 4698|NCTC 2665"         "-" \
    "S_aureus"        "Staphylococcus aureus"       "SA113|ATCC 35556"            "-" \
    "S_capitis"       "Staphylococcus capitis"      "ATCC 27840"                  "-" \
    "S_epidermidis"   "Staphylococcus epidermidis"  "ATCC 12228"                  "-" \
    "S_hominis"       "Staphylococcus hominis"      "ATCC 27844"                  "-" \
    "S_warneri"       "Staphylococcus warneri"      "ATCC 27836"                  "-" \
    "S_mitis"         "Streptococcus mitis"         "ATCC 49456|NCTC 12261"       "-" \
    > "$OUT/cepas.tsv"
fi

# ---------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------
# "A|B" -> "\"A\" OR \"B\"" for Entrez
terms_to_entrez() {
  awk -v t="$1" 'BEGIN{n=split(t,T,"|"); for(i=1;i<=n;i++) printf "%s\"%s\"", (i>1?" OR ":""), T[i]}'
}

# List NCBI assemblies: acc <TAB> level <TAB> name <TAB> species <TAB> infraspecific values
list_assemblies() {  # $1 = Entrez query
  esearch -db assembly -query "$1" </dev/null 2>/dev/null \
    | esummary 2>/dev/null \
    | xtract -pattern DocumentSummary \
        -element AssemblyAccession AssemblyStatus AssemblyName SpeciesName Sub_value
}

# Picks the best assembly whose strain field matches any term.
# Priority: Complete Genome > Chromosome > Scaffold > Contig; GCF_ (RefSeq) before GCA_.
pick_assembly() {  # $1 = species, $2 = terms
  local q; q="\"$1\"[Organism] AND ($(terms_to_entrez "$2")) AND latest[filter]"
  list_assemblies "$q" \
  | awk -F'\t' -v terms="$2" '
      BEGIN{ n=split(terms,T,"|"); for(i=1;i<=n;i++){ t=toupper(T[i]); gsub(/[^A-Z0-9]/,"",t); T[i]=t } }
      { l=""; for(j=2;j<=NF;j++) l=l $j; l=toupper(l); gsub(/[^A-Z0-9]/,"",l);
        for(i=1;i<=n;i++) if(index(l,T[i])>0){ print; next } }' \
  | awk -F'\t' '{ r=($2~/Complete/)?1:($2~/Chromosome/)?2:($2~/Scaffold/)?3:4;
                  g=($1~/^GCF/)?0:1; print r "\t" g "\t" $1 "\t" $2 }' \
  | sort -k1,1n -k2,2n | head -n 1 | cut -f3,4
}

download_genome() {  # $1 = accession, $2 = destination folder
  local acc="$1" d="$2"
  [[ -s "$d/genome.fna" ]] && return 0
  mkdir -p "$d"
  datasets download genome accession "$acc" --include genome \
      --filename "$d/pkg.zip" --no-progressbar >/dev/null 2>&1 || return 1
  unzip -oq "$d/pkg.zip" -d "$d/pkg" || return 1
  cat "$d"/pkg/ncbi_dataset/data/"$acc"/*.fna > "$d/genome.fna" || return 1
  rm -rf "$d/pkg" "$d/pkg.zip"
}

extract_16S() {  # $1 = label, $2 = accession, $3 = genome folder
  local label="$1" acc="$2" d="$3"
  barrnap --kingdom bac --threads "$THREADS" --quiet \
      --outseq "$d/rrna.fa" "$d/genome.fna" > "$d/rrna.gff" 2>/dev/null
  seqkit grep -r -n -p '^16S_rRNA' "$d/rrna.fa" 2>/dev/null \
    | seqkit seq -m "$MINLEN" -M "$MAXLEN" 2>/dev/null \
    | awk -v l="$label" -v a="$acc" \
        '/^>/{ n++; sub(/^>16S_rRNA::/,""); print ">" l "|" a "|copy" n "|" $0; next } { print }' \
    > "$OUT/per_strain/$label.16S.fa"
}

fallback_nuccore() {  # $1 = label, $2 = species, $3 = terms
  local q; q="\"$2\"[Organism] AND ($(terms_to_entrez "$3")) AND 16S ribosomal RNA[Title] AND 1200:1700[SLEN]"
  esearch -db nuccore -query "$q" </dev/null 2>/dev/null | efetch -format fasta 2>/dev/null \
    > "$OUT/fallback_nuccore/$1.nuccore.fa"
  [[ -s "$OUT/fallback_nuccore/$1.nuccore.fa" ]]
}

# ---------------------------------------------------------------------
# Main loop (fd 3 so esearch/efetch don't consume the loop's stdin)
# ---------------------------------------------------------------------
printf 'label\tspecies\tsource\taccession\tlevel\tn_16S_copies\tn_unique_sequences\tlength\n' \
  > "$OUT/resumen.tsv"
: > "$OUT/accesiones_usadas.tsv"

while IFS=$'\t' read -r label species terms manual <&3; do
  [[ -z "${label:-}" || "$label" == \#* ]] && continue
  echo ">>> $label ($species | $terms)"

  acc=""; level=""; fuente="genome"
  if [[ "${manual:--}" != "-" ]]; then
    acc="$manual"; level="manual"
  else
    IFS=$'\t' read -r acc level < <(pick_assembly "$species" "$terms")
  fi

  if [[ -n "$acc" ]] && download_genome "$acc" "$OUT/genomes/$acc"; then
    echo "    genome: $acc ($level)"
    extract_16S "$label" "$acc" "$OUT/genomes/$acc"
  else
    echo "    [!] no genome for this strain; saving species candidates in sin_match/"
    list_assemblies "\"$species\"[Organism] AND \"complete genome\"[filter] AND latest[filter]" \
      | head -n 200 > "$OUT/sin_match/$label.candidatos.tsv"
    acc="-"; level="-"; fuente="nuccore(CHECK)"
    : > "$OUT/per_strain/$label.16S.fa"
    if fallback_nuccore "$label" "$species" "$terms"; then
      awk -v l="$label" '/^>/{ n++; sub(/^>/,""); gsub(/ /,"_"); print ">" l "|nuccore|" n "|" $0; next } { print }' \
        "$OUT/fallback_nuccore/$label.nuccore.fa" > "$OUT/per_strain/$label.16S.fa"
    else
      fuente="NONE"
    fi
  fi

  f="$OUT/per_strain/$label.16S.fa"
  ncop=$(grep -c '>' "$f" || true)
  if (( ncop > 0 )); then
    nuniq=$(seqkit seq -s -w 0 "$f" | sort -u | wc -l)
    lens=$(seqkit fx2tab -n -l "$f" | cut -f2 | sort -n | awk 'NR==1{mn=$1}{mx=$1}END{print mn "-" mx}')
  else
    nuniq=0; lens="-"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$label" "$species" "$fuente" "$acc" "$level" "$ncop" "$nuniq" "$lens" >> "$OUT/resumen.tsv"
  printf '%s\t%s\n' "$label" "$acc" >> "$OUT/accesiones_usadas.tsv"
done 3< "$OUT/cepas.tsv"

# ---------------------------------------------------------------------
# Final FASTA files
# ---------------------------------------------------------------------
cat "$OUT"/per_strain/*.16S.fa > "$OUT/SkinCom_16S_todas_las_copias.fa"
: > "$OUT/SkinCom_16S_1_por_cepa.fa"
for f in "$OUT"/per_strain/*.16S.fa; do
  [[ -s "$f" ]] && seqkit head -n 1 "$f" >> "$OUT/SkinCom_16S_1_por_cepa.fa"
done

echo
echo "=== SUMMARY ==="
column -t -s $'\t' "$OUT/resumen.tsv"
echo
echo "Outputs in $OUT/:"
echo "  SkinCom_16S_1_por_cepa.fa           (1 copy per strain)"
echo "  SkinCom_16S_todas_las_copias.fa     (all copies of the rrn operon)"
echo "  resumen.tsv / accesiones_usadas.tsv (for the methods section)"
echo "  sin_match/, fallback_nuccore/       (only if some strain failed: CHECK manually)"
```

From `SkinCom_16S_1_por_cepa.fa` two derived databases were built: one for `dada2::assignTaxonomy()` (SILVA-like format, down to genus) and a BLAST database (`blastdb/SkinCom_16S`) for species assignment by exact match (see [section 5](#5-taxonomic-assignment-silva--blast)).

---

## 3. Primer trimming (cutadapt)

Primers must be removed before DADA2: their presence interferes with read overlap, the error model, and taxonomic assignment. Amplified region: V3–V4, primers 341F / 805R.

```bash
cd "/mnt/c/AA bioinfo tesis/comunidadmicrobianaepidermidis"

# 1. Sample list
ls raw/*_1.fq.gz | xargs -n 1 basename | sed 's/_1\.fq\.gz$//' > samples.txt
cat samples.txt

# 2. Primer trimming with cutadapt (per sample; -g/-G = forward/reverse primer,
#    -O 15 = minimum overlap, -m 200 = discards pairs shorter than that after
#    trimming, --discard-untrimmed = discards pairs missing either primer)
cutadapt \
  -g ACTCCTACGGGAGGCAGCA \
  -G GGACTACHVGGGTWTCTAAT \
  -O 15 -m 200 --discard-untrimmed -j 0 \
  -o raw/trimmed/${sample}_R1.trimmed.fastq.gz \
  -p raw/trimmed/${sample}_R2.trimmed.fastq.gz \
  raw/${sample}_1.fq.gz raw/${sample}_2.fq.gz \
  > stats/cutadapt_primer_trimming_stats.txt
```

**General principles when adjusting these parameters for another dataset:** check primer orientation, choose amplicon length limits appropriate for the marker used, and inspect the trimming report before continuing — if trimming goes wrong, check primer degeneracy and orientation first. The parameters to adjust are the primer sequences (`-g`/`-G`, must exactly match the pair used) and the `-m`/`-M` limits, which must be adjusted to the expected amplicon length after trimming — don't copy them from another project without checking.

**Verification:** in this dataset, 99.5% of pairs passed the filter (152,433 → 152,433 pairs with both primers detected, 99.6%/99.9% R1/R2), well above the usual >90% threshold.

---

## 4. DADA2: from raw reads to ASVs

```r
library(dada2)

# Define the path where the reads trimmed with Cutadapt are located
path <- "C:/AA Bioinfo tesis/comunidadmicrobianaepidermidis/raw/trimmed"

# Get the paths of the R1 and R2 files
fnFs <- sort(list.files(path, pattern = "_R1.trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.trimmed.fastq.gz", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

# 4. Inspect quality profiles
## These help choose the truncation length; it's one of the most important
## judgment calls in the pipeline, because a poor choice can destroy read
## overlap or leave too many low-quality bases.
if (!dir.exists("stats")) dir.create("stats", recursive = TRUE)

qp_F <- plotQualityProfile(fnFs[1:3])
qp_R <- plotQualityProfile(fnRs[1:3])

ggplot2::ggsave("stats/quality_profile_forward.pdf", qp_F, width = 8, height = 6)
ggplot2::ggsave("stats/quality_profile_reverse.pdf", qp_R, width = 8, height = 6)

qp_F_all <- dada2::plotQualityProfile(fnFs)
qp_R_all <- dada2::plotQualityProfile(fnRs)

ggplot2::ggsave("stats/quality_profile_forward_all.pdf", qp_F_all, width = 16, height = 12)
ggplot2::ggsave("stats/quality_profile_reverse_all.pdf", qp_R_all, width = 16, height = 12)

# 5. Filtering and trimming
## Removes low-quality reads/bases before denoising.
##
## IMPORTANT: truncLen = c(0,0) is intentional, not a placeholder. Since cutadapt
## trims primers tolerating indels, the resulting reads have VARIABLE length
## (225-234 bp in this dataset). A fixed truncLen (e.g. 230) would have
## discarded any read shorter than that value -- in this dataset,
## that would have thrown away ~78% of the reads. maxEE is the real quality filter.
filt_path <- "C:/AA Bioinfo tesis/comunidadmicrobianaepidermidis/raw/trimmed/filtered"
dir.create(filt_path, showWarnings = FALSE)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  trimLeft = c(0, 0),
  truncLen = c(0, 0),       # KEY: 0 so as not to destroy overlap after using Cutadapt
  maxN = 0,
  maxEE = c(2, 2),          # The real quality filter: discards reads with >2 expected errors
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE,      # KEY: FALSE so RStudio doesn't hang on Windows
  matchIDs = TRUE
)

head(out)
write.csv(out, "stats/filterAndTrim_read_tracking.csv")

# 6. Learn the error model
## DADA2 empirically models the sequencing error to distinguish real
## variants from errors.
set.seed(100)

errF <- learnErrors(filtFs, nbases = 1e8, multithread = FALSE, randomize = TRUE)
errR <- learnErrors(filtRs, nbases = 1e8, multithread = FALSE, randomize = TRUE)

ep_F <- plotErrors(errF, nominalQ = TRUE)
ep_R <- plotErrors(errR, nominalQ = TRUE)

ggplot2::ggsave("stats/error_rates_forward.pdf", ep_F, width = 8, height = 6)
ggplot2::ggsave("stats/error_rates_reverse.pdf", ep_R, width = 8, height = 6)

# 7. ASV inference and pair merging
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
head(mergers[[1]])

# 8. Sequence table and chimera removal
seqtab.all <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab.all, method = "consensus", multithread = FALSE)

dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

# Read tracking through the whole pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "sapply_seguimiento_lecturas2.csv", row.names = TRUE)
```

**Result:** 173 final ASVs, average retention of 85.7% (input → non-chimeric) across the 12 samples, with no problematic samples (range 84.8%–87.0%).

---

## 5. Taxonomic assignment (SILVA + BLAST)

Two-step strategy: DADA2's `assignTaxonomy()` down to genus (50% bootstrap), and BLAST against the curated database of the 9 strains to resolve species by exact match.

```r
# 1. Assignment down to genus
taxa <- assignTaxonomy(
  seqtab.nochim,
  refFasta = "refs/SkinCom_16S_dada2_taxonomy.fa",
  multithread = TRUE,
  minBoot = 50
)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

setwd("C:/AA bioinfo tesis/comunidadmicrobianaepidermidis")
library(Biostrings)

# 1) ASVs to FASTA
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV", seq_along(asv_seqs))
writeXStringSet(DNAStringSet(setNames(asv_seqs, asv_ids)), "results/asv.fa")

# 2) BLAST against the curated database (in terminal / system2)
system2("blastn", c("-query", "results/asv.fa", "-db", "blastdb/SkinCom_16S",
                    "-outfmt", shQuote("6 qseqid sseqid pident length qcovs bitscore"),
                    "-max_target_seqs", "50", "-out", "results/asv_vs_db.tsv"))

# 3) Species = best hit (ties are reported together, e.g. "capitis/epidermidis")
hits <- read.delim("results/asv_vs_db.tsv", header = FALSE,
                   col.names = c("asv", "ref", "pident", "len", "qcov", "bitscore"))
hits$especie <- sub("_16S_[0-9]+$", "", hits$ref)

id_min <- 98; qcov_min <- 90; margen_bits <- 6     # 6 bits ~ 1 difference

asignar <- function(d) {
  top <- d[d$bitscore >= max(d$bitscore) - margen_bits, ]
  ok  <- max(top$pident) >= id_min && max(top$qcov) >= qcov_min
  data.frame(asv = d$asv[1],
             Species = if (ok) paste(sort(unique(sub("^[A-Za-z]+_", "", top$especie))), collapse = "/") else NA_character_,
             identidad = max(top$pident), cobertura = max(top$qcov), stringsAsFactors = FALSE)
}
res <- do.call(rbind, lapply(split(hits, hits$asv), asignar))
res <- res[match(asv_ids, res$asv), ]

tax_final <- data.frame(ASV = asv_ids, taxa, Species = res$Species,
                        identidad = res$identidad, cobertura = res$cobertura, row.names = NULL)
write.csv(tax_final, "results/taxonomy_final.csv", row.names = FALSE)

# Consolidation: combine exact match + best BLAST hit into a single
# Species column, with the origin of the assignment recorded.
conf <- !is.na(tax_final$Species) & !is.na(tax_final$Species.1) & tax_final$Species != tax_final$Species.1
sum(conf)  # sanity check: should be 0 or very few
tax_final[conf, c("ASV", "Genus", "Species", "Species.1")]

tax_final$Species_origen <- ifelse(!is.na(tax_final$Species), "exact_match",
                                   ifelse(!is.na(tax_final$Species.1), "best_BLAST_hit", NA_character_))
tax_final$Species   <- ifelse(is.na(tax_final$Species), tax_final$Species.1, tax_final$Species)
tax_final$Species.1 <- NULL

# Genera with only one possible species in the designed community: infer directly
esp_unica <- c(Corynebacterium = "afermentans",
               Micrococcus     = "luteus",
               Cutibacterium   = "acnes")

falta <- is.na(tax_final$Species) & tax_final$Genus %in% names(esp_unica)
tax_final$Species[falta] <- unname(esp_unica[tax_final$Genus[falta]])
tax_final$Species_origen[falta] <- "inferred_from_genus"

table(tax_final$Species, tax_final$Species_origen, useNA = "ifany")
write.csv(tax_final, "results/taxonomy_final.csv", row.names = FALSE)
```

**Result:** 173 ASVs → 7 of the 9 designed species recovered in the final community (*S. hominis* and *M. luteus* absent — expected, with no cross-species mixing artifacts). BLAST identity and coverage of 100% for all species assigned by exact match. A small number of ASVs (e.g. ASV50, ASV72, ASV156) remained unresolved at species level within *Staphylococcus*.

---

## 6. Phylogenetic tree (optional)

Optional step, used only if a UniFrac distance or a tree figure is needed; this work's main beta diversity analysis uses Aitchison/CLR distance (see section 10), which does not require a tree.

```r
library(DECIPHER)
library(phangorn)
library(Biostrings)

seqs <- getSequences(seqtab.nochim)
asv_ids <- paste0("ASV", seq_along(seqs))     # same IDs as in tax_final
names(seqs) <- asv_ids

alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA, processors = NULL)

phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit <- pml(treeNJ, data = phang.align)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(
  fitGTR,
  model = "GTR",
  optInv = TRUE,
  optGamma = TRUE,
  rearrangement = "NNI",
  control = pml.control(trace = 0)
)

tree <- midpoint(fitGTR$tree)                 # root it (needed for UniFrac)

saveRDS(tree, "results/tree.rds")
ape::write.tree(tree, "results/tree.nwk")
```

---

## 7. Building the phyloseq object

```r
# Export ASVs, count table and taxonomy to the formats phyloseq expects
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV", seq_along(asv_seqs))
stopifnot(identical(asv_ids, tax_final$ASV))

write(c(rbind(paste0(">", asv_ids), asv_seqs)), "results/ASVs.fa")

asv_tab <- t(seqtab.nochim)
rownames(asv_tab) <- asv_ids
write.table(asv_tab, "results/ASVs_counts.tsv", sep = "\t", quote = FALSE, col.names = NA)

asv_tax <- tax_final[, setdiff(names(tax_final), "ASV")]
rownames(asv_tax) <- asv_ids
write.table(asv_tax, "results/ASVs_taxonomy.tsv", sep = "\t", quote = FALSE, col.names = NA)

asv_tabtax <- data.frame(asv_tab, asv_tax, check.names = FALSE)
write.table(asv_tabtax, "results/ASVs_counts_taxonomy.tsv", sep = "\t", quote = FALSE, col.names = NA)

# Metadata: from Excel to TSV
library(readxl)
map <- as.data.frame(read_excel("metadata/metadata.xlsx"))
rownames(map) <- map[[1]]
map[[1]] <- NULL
write.table(map, "metadata/metadata.tsv", sep = "\t", quote = FALSE, col.names = NA)

library(phyloseq)
library(ape)
library(Biostrings)

SV <- read.delim("results/ASVs_counts.tsv", row.names = 1, check.names = FALSE)

tax_raw <- read.delim("results/ASVs_taxonomy.tsv", row.names = 1, stringsAsFactors = FALSE)
rangos  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax <- as.matrix(tax_raw[, rangos])

map <- read.delim("metadata/metadata.tsv", row.names = 1, check.names = FALSE)

# Sequencing platform (BMK) IDs -> condition/replicate names
equiv <- c(
  "Bacteria_BMKO260810-AA254-001M0001-01" = "pHD17-1",
  "Bacteria_BMKO260810-AA254-001M0002-01" = "pHD17-2",
  "Bacteria_BMKO260810-AA254-001M0003-01" = "pHD17-3",
  "Bacteria_BMKO260810-AA254-001M0004-01" = "pHD17-4",
  "Bacteria_BMKO260810-AA254-001M0005-01" = "Dbap-1",
  "Bacteria_BMKO260810-AA254-001M0006-01" = "Dbap-2",
  "Bacteria_BMKO260810-AA254-001M0007-01" = "Dbap-3",
  "Bacteria_BMKO260810-AA254-001M0008-01" = "Dbap-4",
  "Bacteria_BMKO260810-AA254-001M0009-01" = "PF-1",
  "Bacteria_BMKO260810-AA254-001M0010-01" = "PF-2",
  "Bacteria_BMKO260810-AA254-001M0011-01" = "PF-3",
  "Bacteria_BMKO260810-AA254-001M0012-01" = "PF-4"
)

stopifnot(setequal(colnames(SV), names(equiv)), !anyDuplicated(equiv))
colnames(SV) <- unname(equiv[colnames(SV)])

stopifnot(setequal(colnames(SV), rownames(map)), setequal(rownames(SV), rownames(tax)))
map <- map[colnames(SV), , drop = FALSE]

ps <- phyloseq(
  otu_table(as.matrix(SV), taxa_are_rows = TRUE),
  sample_data(map),
  tax_table(tax)
)

ps <- merge_phyloseq(ps, readDNAStringSet("results/ASVs.fa"))
if (file.exists("results/tree.nwk")) {
  ps <- merge_phyloseq(ps, phy_tree(read.tree("results/tree.nwk")))
}
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps
```

---

## 8. Quality and prevalence filtering

Technical or very low-information categories that shouldn't be analyzed together with the target community are removed: very low-depth samples and residual ASVs from sequencing noise.

```r
# 0) Depth per sample
sort(sample_sums(ps))

# 1) Samples with very low depth (none falls under this threshold in this dataset)
ps_filtered <- prune_samples(sample_sums(ps) >= 250, ps)
nsamples(ps) - nsamples(ps_filtered)

# 2) Noise: <0.005% of total reads or present in only one sample
umbral <- 0.00005 * sum(sample_sums(ps_filtered))
otu <- as(otu_table(ps_filtered), "matrix")
if (!taxa_are_rows(ps_filtered)) otu <- t(otu)
prevalencia <- rowSums(otu > 0)
stopifnot(identical(names(prevalencia), taxa_names(ps_filtered)))

mantener <- taxa_sums(ps_filtered) >= umbral & prevalencia >= 2
ps_filtered <- prune_taxa(mantener, ps_filtered)

cat("Threshold:", round(umbral), "reads | ASVs:", ntaxa(ps), "->", ntaxa(ps_filtered), "\n")
cat("% reads retained:", round(100 * sum(sample_sums(ps_filtered)) / sum(sample_sums(ps)), 3), "\n")

resumen <- function(x) {
  tt <- as.data.frame(tax_table(x))
  tapply(taxa_sums(x), paste(tt$Genus, tt$Species), sum)
}
a <- resumen(ps); d <- resumen(ps_filtered)
cbind(before = a, after = d[names(a)])
ps_filtered
```

`ps_filtered` is the object used in all downstream analyses (sections 9–13).

---

## 9. Alpha diversity

```r
library(ggpubr)

ps_nonzero <- prune_taxa(taxa_sums(ps_filtered) > 0, ps_filtered)
sample_data(ps_nonzero)$Community <- factor(sample_data(ps_nonzero)$Community, levels = c("pHD17", "Dbap", "PF"))
table(sample_data(ps_nonzero)$Community)      # should come out 4, 4, 4

alfa <- estimate_richness(ps_nonzero, measures = c("Observed", "Shannon", "Simpson"))
rownames(alfa) <- sample_names(ps_nonzero)       # estimate_richness turns hyphens into dots
alfa$Pielou <- alfa$Shannon / log(alfa$Observed)

df <- cbind(as(sample_data(ps_nonzero), "data.frame"), alfa)
df

# With only 4 replicates per group, the test has low power and only detects
# very large differences -- absence of significance does not mean equality.
pares <- combn(levels(df$Community), 2, simplify = FALSE)
grafico <- function(indice) {
  ggboxplot(df, x = "Community", y = indice, fill = "Community", add = "jitter") +
    scale_fill_viridis_d() + theme_classic() + theme(legend.position = "none") +
    stat_compare_means(comparisons = pares, method = "wilcox.test")
}
fig <- ggarrange(grafico("Observed"), grafico("Shannon"), grafico("Simpson"), grafico("Pielou"), ncol = 4)
print(fig)
ggsave("results/alfa_diversidad.pdf", fig, width = 14, height = 4)

kruskal.test(Shannon ~ Community, data = df)
pairwise.wilcox.test(df$Shannon, df$Community, p.adjust.method = "holm")
```

```r
# Half-violin plots (alternative style)
library(see)
library(gghalves)
library(ggplot2)

group_order   <- c("pHD17", "Dbap", "PF")
group_palette <- c(pHD17 = "midnightblue", Dbap = "gold", PF = "orange")

df$Community <- factor(df$Community, levels = group_order)

plot_alfa <- function(indice, etiqueta) {
  ggplot(df, aes(x = Community, y = .data[[indice]])) +
    geom_violinhalf(aes(fill = Community), scale = "width", width = 0.6, flip = TRUE, trim = TRUE) +
    geom_half_boxplot(aes(fill = Community), outlier.size = 0.01, width = 0.3, size = 0.5, alpha = 0.4, side = "r") +
    geom_jitter(width = 0.03, size = 1.5, alpha = 0.8) +
    scale_fill_manual(values = group_palette) +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab(etiqueta) + xlab(NULL)
}

p_shannon  <- plot_alfa("Shannon",  "Shannon diversity")
p_richness <- plot_alfa("Observed", "Observed richness")

ggsave("stats/alpha_shannon.pdf",  p_shannon,  width = 6, height = 6)
ggsave("stats/alpha_richness.pdf", p_richness, width = 6, height = 6)
```

**Result:** no significant differences between pHD17/Dbap/PF (Kruskal-Wallis, Shannon p=0.87) — overall community diversity remains stable; the *bap* effect (see section 13) is species-specific, not a general shift in diversity.

---

## 10. Beta diversity (Aitchison/CLR)

Aitchison/CLR distance instead of Bray-Curtis/UniFrac: more appropriate for 16S compositional data (Gloor et al. 2017), and does not require a phylogenetic tree.

```r
library(vegan)
library(zCompositions)

X <- as(otu_table(ps_filtered), "matrix")
if (taxa_are_rows(ps_filtered)) X <- t(X)
dim(X)
sum(X == 0)

d.czm <- cmultRepl(X, method = "CZM", label = 0)   # zero replacement
E.clr <- t(apply(d.czm, 1, function(x) log(x) - mean(log(x))))   # CLR per sample

d.pcx    <- prcomp(E.clr)
dist.clr <- dist(E.clr)                            # Aitchison distance
summary(d.pcx)$importance[, 1:4]

# PCA on CLR (with pseudocount 1, within phyloseq)
group_palette <- c(pHD17 = "midnightblue", Dbap = "gold", PF = "orange")

ps_clr <- transform_sample_counts(ps_filtered, function(x) { x <- x + 1; log(x) - mean(log(x)) })
ord    <- ordinate(ps_clr, "RDA", "euclidean")

pca_plot <- plot_ordination(ps_clr, ord, color = "Community") +
  theme_classic() +
  geom_point(aes(color = Community), alpha = 0.8, size = 3) +
  scale_colour_manual(values = group_palette) +
  stat_ellipse(lwd = 1)

print(pca_plot)
ggsave("results/beta_pca_aitchison.pdf", pca_plot, width = 6, height = 6)

# PERMANOVA on the same CLR (Aitchison) distances
X <- as(otu_table(ps_clr), "matrix")
if (taxa_are_rows(ps_clr)) X <- t(X)
clr_dist <- dist(X, method = "euclidean")

meta_df <- data.frame(sample_data(ps_clr))
meta_df <- meta_df[rownames(X), , drop = FALSE]
meta_df$Community <- as.factor(meta_df$Community)

set.seed(123)
adonis_res <- adonis2(clr_dist ~ Community, data = meta_df, permutations = 999)
print(adonis_res)

# Dispersion homogeneity across groups (PERMANOVA assumption)
anova(betadisper(clr_dist, meta_df$Community))

# ANOSIM: non-parametric confirmation based on distance ranks
set.seed(123)
ano.group <- anosim(clr_dist, meta_df$Community, permutations = 999)
ano.group

pdf("results/beta_anosim.pdf", width = 6, height = 6)
plot(ano.group)
dev.off()
```

**Result:** R = 0.007, P = 0.398 (ANOSIM) — the three communities do not separate in ordination space. With n=4/group, power to separate centroids is limited; consistent with the rest of the analysis (the *bap* effect is species-specific, not a global compositional shift).

---

## 11. Rarefaction curves

```r
library(vegan)
library(ggplot2)

X <- as(otu_table(ps_filtered), "matrix")
if (taxa_are_rows(ps_filtered)) X <- t(X)
meta <- data.frame(sample_data(ps_filtered))

set.seed(123)
n_rep <- 10
curvas <- do.call(rbind, lapply(rownames(X), function(s) {
  x <- X[s, ]; N <- sum(x)
  prof <- unique(c(100, 250, 500, seq(1000, N, by = 2000), N))
  data.frame(muestra = s, lecturas = prof,
             shannon = sapply(prof, function(d)
               mean(replicate(n_rep, suppressWarnings(
                 vegan::diversity(vegan::rrarefy(x, d), index = "shannon"))))))
}))
curvas$Community <- factor(meta[curvas$muestra, "Community"], levels = c("pHD17", "Dbap", "PF"))

group_palette <- c(pHD17 = "midnightblue", Dbap = "gold", PF = "orange")
p <- ggplot(curvas, aes(lecturas, shannon, group = muestra, color = Community)) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(labels = scales::comma) +
  scale_color_manual(values = group_palette) +
  theme_classic() + xlab("Reads") + ylab("Shannon")

print(p)
ggsave("stats/shannon_rarefaction.pdf", p, width = 7, height = 5)
```

---

## 12. Composition: bar plots, bubble plot, heatmap, Venn, UpSet

### 12.1 Relative abundance bar plots (by sample and by community)

```r
library(dplyr)
library(ggplot2)

nivel <- "Species"

ps_rel <- ps_filtered %>%
  tax_glom(taxrank = nivel, NArm = FALSE) %>%
  transform_sample_counts(function(x) x / sum(x))

# By sample
Sample_tax <- psmelt(ps_rel) %>% arrange(.data[[nivel]])
Sample_tax$Sample <- factor(Sample_tax$Sample, levels = sample_names(ps_filtered))

p1 <- ggplot(Sample_tax, aes(x = Sample, y = Abundance, fill = .data[[nivel]])) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Community, scales = "free_x", space = "free_x") +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Relative abundance") + ggtitle(paste(nivel, "composition by sample"))
ggsave("stats/barras_muestra.pdf", p1, width = 10, height = 5)

# By community (average of replicates)
type_tax <- Sample_tax %>%
  group_by(Community, .data[[nivel]]) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

p2 <- ggplot(type_tax, aes(x = Community, y = Abundance, fill = .data[[nivel]])) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_classic() + ylab("Relative abundance") + xlab(NULL) +
  ggtitle(paste(nivel, "composition by community (mean of replicates)"))
ggsave("stats/barras_comunidad.pdf", p2, width = 6, height = 5)

# Percentage table (the bars don't show minor species well)
tt <- as.data.frame(type_tax); tt$taxon <- tt[[nivel]]
round(100 * xtabs(Abundance ~ taxon + Community, data = tt), 2)
```

### 12.2 Bubble plot (mean ± SD by species and community)

```r
# ================================================================
# Bubble plot of relative abundance by species, with SD across
# replicates (object ps_filtered, Community column with levels
# pHD17 / Dbap / PF)
# ================================================================
library(phyloseq)
library(dplyr)
library(data.table)
library(ggplot2)

nivel <- "Species"

# --- 1) SD across replicates, by species and community ----------
# Calculated BEFORE merging samples: each replicate is normalized to
# % on its own, and then SD is calculated within each Community.
sd1 <- ps_filtered %>%
  tax_glom(taxrank = nivel, NArm = FALSE) %>%
  transform_sample_counts(function(x) x / sum(x) * 100) %>%
  psmelt() %>%
  group_by(OTU, Community) %>%
  summarize(SD = sd(Abundance, na.rm = TRUE), .groups = "drop")

# --- 2) Mean by community --------------------------------------
bub1 <- merge_samples(ps_filtered, "Community")

type_taxa1 <- bub1 %>%
  tax_glom(taxrank = nivel, NArm = FALSE) %>%
  transform_sample_counts(function(x) x / sum(x) * 100) %>%
  psmelt() %>%
  mutate(Community = as.character(Sample)) %>%   # merge_samples scrambles metadata; rebuild Community
  arrange(OTU)

dat1 <- type_taxa1 %>%
  left_join(sd1, by = c("OTU", "Community")) %>%
  as.data.table()

# --- 3) Rare species are NOT collapsed -----------------------------
# Unlike a noisy environmental community with hundreds of taxa,
# this is a DEFINED synthetic community of 9 species: if a species
# drops sharply in one condition (e.g. the Dbap background), that's
# exactly what should be shown, not hidden. Only the ~0.01% of reads
# with no species assignment is collapsed.
dat1[, Species_plot := as.character(get(nivel))]
dat1[is.na(Species_plot) | Species_plot == "NA", Species_plot := "Unassigned (background noise, ~0.01%)"]

# Order of the Y axis and legend: by descending mean abundance.
# Assigning the SAME factor levels before plotting keeps the legend's
# color order in sync with the Y-axis bubbles.
orden_especies <- dat1[, .(m = mean(Abundance, na.rm = TRUE)), by = Species_plot][order(-m), Species_plot]
dat1[, Species_plot := factor(Species_plot, levels = orden_especies)]
sample_order <- c("pHD17", "Dbap", "PF")

# --- 4) Plot -----------------------------------------------------
Fig_bubble <- ggplot(
  dat1[Abundance > 0],
  aes(x = factor(Community, levels = sample_order), y = Species_plot)
) +
  geom_point(aes(size = Abundance, fill = Species_plot, color = Species_plot),
             alpha = 0.5, shape = 21) +
  geom_point(aes(size = Abundance + SD, fill = Species_plot, color = Species_plot),
             alpha = 0.2, shape = 21) +
  geom_point(aes(size = pmax(Abundance - SD, 0), fill = Species_plot, color = Species_plot),
             alpha = 1, shape = 21) +
  scale_size_continuous(limits = c(0.000001, 100), range = c(0.1, 12),
                        breaks = c(1, 5, 10, 25, 50)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_bw() +
  labs(x = NULL, y = NULL, size = "Relative abundance (%)",
       fill = "Species", color = "Species",
       title = "Species composition by community 
       (mean \u00b1 SD across replicates)")

print(Fig_bubble)
ggsave("results/bubble_especies.pdf", Fig_bubble, width = 4, height = 7)
```

### 12.3 Heatmap (ampvis2)

```r
# ================================================================
# Heatmap of the most abundant species per community, with ampvis2
# ================================================================
library(ampvis2)
library(phyloseq)

ampvis2_obj <- amp_load(ps_filtered)
n_species <- length(unique(tax_table(ps_filtered)[, "Species"]))

xx <- amp_heatmap(
  ampvis2_obj,
  group_by = "Community",
  tax_aggregate = "Species",   # aggregates to species, not OTU
  tax_add = "Genus",           # adds the genus next to the species name on the Y axis
  tax_show = n_species,        # shows all defined species, none hidden
  color_vector = c("snow", "steelblue"),
  plot_values = TRUE
) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Most abundant species per community")

print(xx)
ggplot2::ggsave("results/heatmap_especies.pdf", xx, width = 5, height = 6)
```

### 12.4 Venn diagram (ampvis2)

```r
# ================================================================
# Venn diagram of species shared/unique across communities
# ================================================================
library(ampvis2)
library(phyloseq)

# amp_venn has no aggregation parameter like amp_heatmap's tax_aggregate,
# so it's aggregated to species first (tax_glom), same as in bubble/heatmap.
ps_species <- tax_glom(ps_filtered, taxrank = "Species", NArm = FALSE)
ampvis2_species <- amp_load(ps_species)

venn <- amp_venn(
  ampvis2_species,
  group = "Community",
  cut_a = 0.01,   # minimum mean relative abundance (%) to count a species as "present"
  cut_f = 75,     # minimum % of the group's replicates in which it must appear (75% = 3 of 4)
  detailed_output = TRUE
)

venn$plot
```

**Result:** the 6 species with abundance >0.01% are shared by all 3 conditions (100% core) — no species exclusive to a single condition.

### 12.5 UpSet plot

```r
# ================================================================
# UpSet plot of species shared/unique across communities
# ================================================================
library(UpSetR)
library(phyloseq)

ps_species <- tax_glom(ps_filtered, taxrank = "Species", NArm = FALSE)

abund_table <- as.data.frame(otu_table(ps_species))
if (!taxa_are_rows(ps_species)) abund_table <- as.data.frame(t(abund_table))
metadata <- data.frame(sample_data(ps_species))

rownames(abund_table) <- as.character(tax_table(ps_species)[rownames(abund_table), "Species"])
rownames(abund_table)[is.na(rownames(abund_table))] <- "Unassigned"

groups <- c("pHD17", "Dbap", "PF")

species_list <- setNames(
  lapply(groups, function(g) {
    cols <- metadata[["Community"]] == g
    rownames(abund_table)[rowSums(abund_table[, cols, drop = FALSE]) > 0]
  }),
  groups
)

upset(fromList(species_list),
      nsets = length(species_list),
      order.by = "freq",
      decreasing = TRUE,
      point.size = 4,
      line.size = 1.2,
      mainbar.y.label = "Number of shared species",
      sets.x.label = "Total species per community",
      text.scale = c(1.5, 1.2, 1.2, 1.2, 1.5, 1.1),
      mb.ratio = c(0.6, 0.4))
```

With only 3 groups and a 100% core (see Venn), the UpSet plot adds no information beyond the Venn diagram — it's included here for completeness but was dropped from the final results presentation.

---

## 13. Differential abundance: LEfSe, ANOVA, ANCOM, ANCOM-BC2

Four independent methods, with different statistical logic, to assess whether the functional *bap* gene confers an advantage to any species in the community. All four are run on `ps_filtered` aggregated to species level (`tax_glom`).

### 13.1 LEfSe

> **Reproducibility note:** `microbiomeMarker`'s `run_lefse()` (v1.13.2, Bioconductor devel) has a real, confirmed bug in its own source code — see [section 14](#14-known-package-issues-and-applied-fixes). The code below uses a patched copy, `run_lefse_fixed()`.

```r
# ================================================================
# LEfSe across the three communities, with microbiomeMarker
# ================================================================
library(microbiomeMarker)
library(phyloseq)

# PATCH: see section 14 for the full explanation of the bug and the fix.
run_lefse_fixed <- function(ps, group, subgroup = NULL, taxa_rank = "all",
                             transform = c("identity", "log10", "log10p"),
                             norm = "CPM", norm_para = list(), kw_cutoff = 0.05,
                             lda_cutoff = 2, bootstrap_n = 30, bootstrap_fraction = 2/3,
                             wilcoxon_cutoff = 0.05, multigrp_strat = FALSE,
                             strict = c("0", "1", "2"), sample_min = 10,
                             only_same_subgrp = FALSE, curv = FALSE) {
    if (!inherits(ps, "phyloseq")) {
        stop("`ps` must be phyloseq object", call. = FALSE)
    }
    ps <- check_rank_names(ps)
    ps <- check_taxa_rank(ps, taxa_rank)
    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    strict <- match.arg(strict, c("0", "1", "2"))
    strict <- as.numeric(strict)
    summarized <- check_tax_summarize(ps)
    if (summarized && norm != "CPM") {
        stop("`norm` must be a 'CPM' or 'none' while `ps` has been summarized",
            call. = FALSE)
    }
    ps <- preprocess_ps(ps)
    ps <- transform_abundances(ps, transform = transform)
    norm_para <- c(norm_para, method = norm, object = list(ps))
    ps_normed <- do.call(normalize, norm_para)
    sample_meta <- sample_data(ps_normed)
    grp_info <- lefse_format_grp(sample_meta, group, subgroup = subgroup)
    grp <- grp_info$group
    subgrps <- grp_info$subgroup
    grp_hie <- grp_info$group_hie
    ps_summarized <- pre_ps_taxa_rank(ps_normed, taxa_rank)
    otus <- abundances(ps_summarized, norm = TRUE)
    otus_test <- as.data.frame(t(otus), stringsAsFactors = FALSE)
    feature <- tax_table(ps_summarized)@.Data[, 1]
    names(otus_test) <- feature
    tax <- matrix(feature) %>% tax_table()
    row.names(tax) <- row.names(otus)
    kw_p <- purrr::map_dbl(otus_test, ~kruskal.test(.x, grp)$p.value)
    na_ind <- is.na(kw_p)
    if (sum(na_ind) >= 1) {
        otus_test <- otus_test[!na_ind]
        kw_p <- kw_p[!na_ind]
    }
    sig_ind <- kw_p <= kw_cutoff
    sig_otus <- otus_test[, sig_ind, drop = FALSE]   # <-- THE FIX (was: otus_test[, sig_ind])
    features_nms <- names(sig_otus)
    wilcoxon_p <- purrr::map2_lgl(sig_otus, features_nms, ~test_rep_wilcoxon(subgrps,
        grp_hie, .x, .y, wilcoxon_cutoff = wilcoxon_cutoff, multicls_strat = multigrp_strat,
        strict = strict, sample_min = sample_min, only_same_subcls = only_same_subgrp,
        curv = curv))
    sig_otus <- sig_otus[, wilcoxon_p, drop = FALSE]
    if (ncol(sig_otus) == 0) {
        warning("No marker was identified", call. = FALSE)
        mm <- microbiomeMarker(marker_table = NULL, norm_method = get_norm_method(norm),
            diff_method = "lefse", otu_table = otu_table(otus,
                taxa_are_rows = TRUE), sam_data = sample_data(ps_normed),
            tax_table = tax)
        return(mm)
    }
    otus_enriched_group <- get_feature_enrich_group(grp, sig_otus)
    ldas <- bootstap_lda(sig_otus, boot_n = bootstrap_n, class = grp,
        sample_fract = bootstrap_fraction)
    lefse_res <- data.frame(feature = names(sig_otus), enrich_group = otus_enriched_group$group,
        ef_lda = ldas, pvalue = kw_p[sig_ind][wilcoxon_p], stringsAsFactors = FALSE)
    lefse_sig <- filter(lefse_res, .data$ef_lda >= lda_cutoff) %>%
        arrange(.data$enrich_group, desc(.data$ef_lda))
    lefse_out <- return_marker(lefse_sig, lefse_res)
    lefse_out$padj <- lefse_out$pvalue
    mm <- microbiomeMarker(marker_table = lefse_out, norm_method = get_norm_method(norm),
        diff_method = "lefse", otu_table = otu_table(otus, taxa_are_rows = TRUE),
        sam_data = sample_data(ps_normed), tax_table = tax)
    mm
}
environment(run_lefse_fixed) <- asNamespace("microbiomeMarker")

ps_species <- tax_glom(ps_filtered, taxrank = "Species", NArm = FALSE)

mm <- run_lefse_fixed(ps_species, group = "Community",
                       norm = "CPM", multigrp_strat = TRUE,
                       sample_min = 4)

head(marker_table(mm))
n_markers <- nrow(marker_table(mm))

p_lefse <- plot_ef_bar(mm, label_level = 1) +
  scale_fill_manual(values = group_palette) +
  scale_y_discrete(labels = function(x) sub("^s__", "", x)) +
  labs(title = "LEfSe markers by community",
       subtitle = paste0(n_markers, " marker(s) -- Kruskal-Wallis p <= 0.05, LDA >= 2"),
       x = "LDA score (log10)", y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40")
  )

print(p_lefse)
ggsave("results/lefse_markers.pdf", p_lefse, width = 15, height = max(2.5, 0.7 * n_markers + 1.3))
```

**Result:** 1 marker — *S. epidermidis*, enriched in pHD17 (Kruskal-Wallis p≤0.05, LDA≥2).

### 13.2 ANOVA (parametric, on CLR)

> **Reproducibility note:** `microbiomeMarker`'s `run_test_multiple_groups()` was dropped here — its `transform` parameter doesn't accept `"clr"`, and it routes through the same internal `preprocess_ps()` that breaks with already CLR-transformed data (see section 14). Used instead: a straightforward, transparent per-species one-way ANOVA with `stats::aov()`.

```r
library(phyloseq)
library(microbiome)

ps_species <- tax_glom(ps_filtered, taxrank = "Species", NArm = FALSE)
ps_clr <- microbiome::transform(ps_species, "clr", shift = 1)

otu <- as(otu_table(ps_clr), "matrix")
if (taxa_are_rows(ps_clr)) otu <- t(otu)

tax <- as(tax_table(ps_species), "matrix")
colnames(otu) <- tax[colnames(otu), "Species"]

meta <- as(sample_data(ps_clr), "data.frame")
group <- factor(meta$Community, levels = c("pHD17", "Dbap", "PF"))
stopifnot(identical(rownames(otu), rownames(meta)))

anova_results <- do.call(rbind, lapply(colnames(otu), function(sp) {
  fit <- aov(otu[, sp] ~ group)
  tab <- summary(fit)[[1]]
  data.frame(
    Species  = sp,
    F_value  = tab["group", "F value"],
    p_value  = tab["group", "Pr(>F)"],
    row.names = NULL
  )
}))

anova_results$p_adj_fdr <- p.adjust(anova_results$p_value, method = "fdr")
anova_results <- anova_results[order(anova_results$p_value), ]

markers <- anova_results
markers
write.csv(markers, "results/markers_anova.csv", row.names = FALSE)
```

**Result:** no species survives FDR correction (q<0.05); *S. epidermidis* comes closest (F=5.85, p=0.024, q=0.165).

### 13.3 ANCOM

```r
# ================================================================
# ANCOM (original, Mandal et al. 2015)
# run_ancom() does its OWN internal compositional handling (confirmed
# via args(run_ancom): transform only accepts identity/log10/log10p,
# no "clr", and has its own norm="TSS" by default) -- it's given
# raw ps_species, not the CLR-transformed object.
# ================================================================
library(microbiomeMarker)
library(phyloseq)

ps_species <- tax_glom(ps_filtered, taxrank = "Species", NArm = FALSE)

mm_ancom <- run_ancom(
  ps_species,
  group = "Community",
  taxa_rank = "Species",
  norm = "TSS",
  p_adjust = "none",
  W_cutoff = 0.75
)

ancom <- marker_table(mm_ancom)
ancom
ancom_df <- as(ancom, "data.frame")   # needed: write.csv() directly on the S4 object fails (see section 14)
write.csv(ancom_df, "results/markers_ANCOM.csv", row.names = FALSE)
```

**Result:** no marker identified at the standard threshold (W_cutoff=0.75). With only 7 species, ANCOM requires a species to be significant against 75% of the others in pairwise comparison — a very high bar for so few taxa. Consistent with the ANOVA-FDR result.

### 13.4 ANCOM-BC2 (combined model, 3 groups)

```r
# ================================================================
# ANCOM-BC2 across the three communities
# No real covariate available in metadata.tsv (only Community, an
# id, and DNA ng/uL -- a technical extraction value, not included).
# Reference = Dbap: same plasmid as pHD17, mutated bap -- isolates
# the bap effect (pHD17 vs Dbap) from the plasmid-burden effect
# (PF vs Dbap).
# ================================================================
library(ANCOMBC)
library(phyloseq)

sample_data(ps_filtered)$Community <- factor(
  sample_data(ps_filtered)$Community,
  levels = c("Dbap", "PF", "pHD17")
)

output_combined <- ancombc2(
  data = ps_filtered,
  assay_name = "counts",
  tax_level = "Species",
  fix_formula = "Community",
  rand_formula = NULL,
  p_adj_method = "fdr",
  pseudo = 0,
  pseudo_sens = TRUE,
  prv_cut = 0.01,       # deliberately lax: only 7-9 designed species, don't accidentally filter any out
  lib_cut = 0,
  s0_perc = 0.05,
  group = "Community",
  struc_zero = TRUE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  trend = FALSE   # no natural order between pHD17/Dbap/PF; trend=TRUE requires contrast matrices that aren't needed
)

res <- output_combined$res
res
write.csv(res, "results/markers_ANCOMBC2.csv", row.names = FALSE)
```

> **Important reproducibility note:** with `pairwise=TRUE` + `dunnet=TRUE` together, `output_combined$res_pair` has broken `p_`/`q_` columns (see section 14) — for pairwise comparisons, use the version in section 13.5 instead. `res` (without pairwise) IS reliable for the two reference contrasts (pHD17 vs Dbap, PF vs Dbap).

### 13.5 ANCOM-BC2 (3 clean pairwise comparisons)

```r
# ================================================================
# ANCOM-BC2, as 3 independent two-group comparisons instead of the
# combined pairwise+dunnet run in section 13.4.
#
# Why: with global/pairwise/dunnet=TRUE all together, res_pair's p/q
# columns come out broken (p=1 for everything, even a case with
# W=3.4 that should give p=0.009 -- a real, confirmed bug, not ours).
# Subsetting to two groups and running plain ancombc2 (pairwise=FALSE
# etc.) is the standard workaround for multi-level factors and
# completely avoids that broken code path.
# ================================================================
library(ANCOMBC)
library(phyloseq)
library(dplyr)

run_pair <- function(ps, level_ref, level_other, contrast_label) {
  # prune_samples(), not subset_samples(): subset_samples() uses NSE
  # internally and can't find level_ref/level_other when called
  # from inside a function (a real phyloseq limitation, not our bug)
  keep <- sample_data(ps)$Community %in% c(level_ref, level_other)
  ps_pair <- prune_samples(keep, ps)
  sample_data(ps_pair)$Community <- droplevels(factor(
    sample_data(ps_pair)$Community, levels = c(level_ref, level_other)
  ))

  fit <- ancombc2(
    data = ps_pair,
    assay_name = "counts",
    tax_level = "Species",
    fix_formula = "Community",
    rand_formula = NULL,
    p_adj_method = "fdr",
    pseudo = 0,
    pseudo_sens = TRUE,
    prv_cut = 0.01,
    lib_cut = 0,
    s0_perc = 0.05,
    group = "Community",
    struc_zero = TRUE,
    neg_lb = FALSE,
    alpha = 0.05,
    n_cl = 2,
    verbose = TRUE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE
  )

  res_col <- paste0("Community", level_other)
  fit$res %>%
    dplyr::transmute(
      taxon,
      lfc = .data[[paste0("lfc_", res_col)]],
      se  = .data[[paste0("se_", res_col)]],
      p   = .data[[paste0("p_", res_col)]],
      q   = .data[[paste0("q_", res_col)]],
      contrast = contrast_label
    )
}

set.seed(321)
df_phd   <- run_pair(ps_filtered, "Dbap", "pHD17", "pHD17 vs Dbap")
df_pf    <- run_pair(ps_filtered, "Dbap", "PF",    "PF vs Dbap")
df_phdpf <- run_pair(ps_filtered, "PF",   "pHD17", "pHD17 vs PF")

df_res <- dplyr::bind_rows(df_phd, df_pf, df_phdpf) %>%
  dplyr::mutate(
    contrast = factor(contrast, levels = c("pHD17 vs Dbap", "PF vs Dbap", "pHD17 vs PF")),
    sig = dplyr::case_when(
      q < 0.05 ~ "q < 0.05 (FDR)",
      q < 0.10 ~ "q < 0.10 (suggestive)",
      TRUE     ~ "ns"
    )
  )

print(df_res)
write.csv(df_res, "results/markers_ANCOMBC2_pairwise.csv", row.names = FALSE)
```

**Result:** pHD17 vs Dbap: lfc=1.24, q=0.095 (the closest to significance, *S. epidermidis*). PF vs Dbap and pHD17 vs PF: no significant or suggestive markers.

### 13.6 ANCOM-BC2 visualization (bar plot + volcano, 3 contrasts)

```r
# ================================================================
# ANCOM-BC2 visualizations: all 3 pairwise contrasts
# Requires section 13.5 to have been run first -- uses its `df_res`.
# ================================================================
library(dplyr)
library(ggplot2)

orden_taxa <- df_res %>%
  dplyr::filter(contrast == "pHD17 vs Dbap") %>%
  dplyr::arrange(lfc) %>%
  dplyr::pull(taxon)
df_res$taxon <- factor(df_res$taxon, levels = orden_taxa)

# ---- 1. Log fold change bar plot, 3 contrasts ----
p_bar <- ggplot(df_res, aes(x = taxon, y = lfc, fill = lfc)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.2) +
  geom_text(
    data = dplyr::filter(df_res, sig != "ns"),
    aes(y = lfc + sign(lfc) * (se + 0.25), label = "*"),
    size = 6, vjust = 0.75
  ) +
  coord_flip() +
  facet_wrap(~contrast, nrow = 1) +
  scale_fill_viridis_c(option = "plasma", name = "LFC") +
  labs(
    title = "ANCOM-BC2: log fold change across all 3 pairwise contrasts",
    subtitle = "* = q < 0.10 (FDR); none reach q < 0.05",
    x = NULL, y = "Log fold change"
  ) +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(face = "bold"), legend.position = "right")

print(p_bar)
ggsave("results/ancombc2_lfc_barplot.pdf", p_bar, width = 11, height = 4.5)

# ---- 2. Volcano plot, 3 contrasts ----
p_volcano <- ggplot(df_res, aes(x = lfc, y = -log10(q))) +
  geom_point(aes(color = sig), size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.10), linetype = "dotted", color = "orange") +
  geom_text(
    data = dplyr::filter(df_res, sig != "ns"),
    aes(label = taxon), hjust = 0.5, vjust = 1.6, size = 3.2
  ) +
  facet_wrap(~contrast, nrow = 1) +
  scale_color_manual(values = c(
    "q < 0.05 (FDR)"        = "red",
    "q < 0.10 (suggestive)" = "orange",
    "ns"                    = "grey60"
  )) +
  labs(
    title = "ANCOM-BC2 volcano: effect size vs statistical support, all 3 contrasts",
    x = "Log fold change",
    y = expression(-log[10](q)),
    color = NULL
  ) +
  theme_minimal(base_size = 13)

print(p_volcano)
ggsave("results/ancombc2_volcano.pdf", p_volcano, width = 11, height = 4.5)
```

### 13.7 Summary of the 4 methods

| Method | Multiple-testing correction? | Result for *S. epidermidis* |
|---|---|---|
| LEfSe | No (KW + Wilcoxon + LDA) | Only marker, enriched in pHD17 |
| ANOVA (CLR) | Yes (FDR) | p=0.024, q=0.165 — closest, not significant |
| ANCOM | Implicit (W_cutoff) | No marker at the standard threshold |
| ANCOM-BC2 (pairwise) | Yes (FDR) | pHD17 vs Dbap: q=0.095 — closest |

Additionally validated independently against the commercial sequencing report (Metastats/BMK): raw p<0.05 in all 3 pairwise comparisons for *S. epidermidis*, but q>0.05 in all 3 — the exact same pattern, confirming neither analysis contradicts the other, they simply differ in whether raw p or corrected q is highlighted as the "significant" criterion.

**Conclusion:** no method reaches formal significance after correction for multiple comparisons (n=4/group limits power), but the directional convergence across methods with different statistical logic consistently points to *S. epidermidis* as the species most affected by *bap* functionality.

---

## 14. Known package issues and applied fixes

Documented here because these are real bugs in the installed versions (not usage errors), reproducible, and relevant to anyone reusing this code with the same package versions.

### 14.1 `microbiomeMarker::run_lefse()` — bare vector when selecting a single column

**Package/version:** `microbiomeMarker` 1.13.2 (Bioconductor devel).

**Symptom:** `Error: Can't recycle .x (size N) to match .y (size 0)` inside `purrr::map2_lgl()`, intermittently depending on how many taxa pass the Kruskal-Wallis filter.

**Root cause:** in `run_lefse()`'s source code, the line `sig_otus <- otus_test[, sig_ind]` uses R's default `data.frame` subsetting (`drop = TRUE`), which converts the result to a flat vector when exactly **one** column passes the filter — breaking the following call to `purrr::map2_lgl()`, which expects a list/data.frame.

**Fix:** an exact copy of the function with `drop = FALSE` added on that line (`run_lefse_fixed()`, section 13.1), with `environment()` reassigned to the package's namespace so it can still see its non-exported internal functions (`preprocess_ps`, `test_rep_wilcoxon`, etc.).

### 14.2 Internal `preprocess_ps()` — misreads already CLR-transformed data

**Symptom:** `The library size of sample(s): X is/are zero` when passing an object already CLR-transformed (negative values) to `run_lefse()` or `run_test_multiple_groups()`.

**Root cause:** the internal helper `preprocess_ps()`, shared by several of the package's `run_*` functions, computes "library size" as the sum of values per sample — an operation that assumes non-negative raw counts. With CLR data (mean ~0 per sample, half the values negative), that sum comes out near zero, and the function interprets it as zero sequencing depth.

**Fix:** never pre-transform to CLR before passing data to these functions — let their own `norm`/`transform` parameter do the normalization internally (`run_lefse(norm="CPM")`, `run_ancom(norm="TSS")`). For the ANOVA, avoid the package's function entirely and use `stats::aov()` directly on the CLR matrix (section 13.2).

### 14.3 `ancombc2()` — degenerate `p`/`q` in `res_pair` with `pairwise` + `dunnet` combined

**Package/version:** `ANCOMBC` (Bioconductor release, not devel — unlike the two above).

**Symptom:** with `pairwise=TRUE` and `dunnet=TRUE` simultaneously, every `p_`/`q_` column in `output_combined$res_pair` comes out exactly `1`, regardless of the actual `W` statistic value.

**Confirmed root cause:** comparing `res` (simple model, no pairwise) row-by-row against `res_pair` for the same nominal coefficient, `lfc`, `se` and `W` are **identical** between the two objects (e.g. W=3.397898 in both), but `p` goes from 0.0094 (correct) to 1 (impossible for that W). The underlying model is the same; something in the joint p-value computation of the `pairwise+dunnet` branch with small n (n=4/group here) returns a degenerate value.

**Fix:** don't use `res_pair` from a combined run. Instead, for each pair of groups, run `ancombc2()` independently and simply (`pairwise=FALSE`, `dunnet=FALSE`, `global=FALSE`) on a subset of just those two groups (`prune_samples()`, section 13.5) — exactly as reliable as the reference run that did work correctly.

### 14.4 `phyloseq::subset_samples()` — non-standard evaluation broken inside a function

**Symptom:** `error evaluating argument 'table' ... object 'level_ref' not found`, when calling `subset_samples(ps, Community %in% c(level_ref, level_other))` inside a user-defined function.

**Root cause:** `subset_samples()` uses non-standard evaluation (NSE) internally and does not correctly resolve variables defined in the calling function's environment — a known phyloseq limitation, not specific to this analysis.

**Fix:** use `prune_samples()` with an explicitly computed logical vector (`sample_data(ps)$Community %in% c(...)`) instead of `subset_samples()`, whenever the subsetting happens inside a function.

### 14.5 `write.csv()` on an S4 `marker_table` object with rows

**Symptom:** `Error in validObject(.Object): invalid class "marker_table" object: marker table must have non-zero dimensions`, when calling `write.csv()` directly on the result of `marker_table()` when there ARE markers (with zero markers, `marker_table()` returns `NULL` and there's no error — different behavior).

**Root cause:** `marker_table` is an S4 class; `write.csv()` internally tries to coerce it to `data.frame`, and that coercion path re-triggers the `marker_table` S4 class validation, which fails at that point in `write.table()`'s pipeline.

**Fix:** explicitly convert with `as(x, "data.frame")` before `write.csv()` (section 13.3) — completely avoids the problematic internal coercion path.

---

## 15. Reproducibility

Still to be included in the final version published alongside the manuscript:

```r
sessionInfo()
# or, if the project uses renv:
# renv::snapshot()
```

**Summary of relevant versions used in this analysis:**
- `microbiomeMarker` 1.13.2 (Bioconductor devel) — see bugs documented in 14.1–14.2
- `ANCOMBC` (Bioconductor release) — see bug documented in 14.3
- `phyloseq`, `dada2`, `vegan`, `ampvis2`, `UpSetR`, `zCompositions`, `compositions`, `microbiome` — exact version pending `sessionInfo()`

**Input data:** 12 samples (4 replicates × 3 conditions: pHD17, Dbap, PF), 16S V3–V4 sequencing (BMK, project BMKO260810-AA254-01-01). Reference database: 9 SkinCom strains (Lekbua et al. 2024).
