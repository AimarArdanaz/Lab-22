# Pipeline de análisis 16S rRNA — Comunidad sintética 

Código completo y reproducible del análisis de la comunidad microbiana sintética (9 especies, basada en Lekbua et al. 2024, *Cell Rep Methods* 4:100832), comparando tres condiciones experimentales: **pHD17** (cepa con plásmido, gen *bap* funcional), **Dbap** (cepa con plásmido, *bap* mutado) y **PF** (cepa sin plásmido). 4 réplicas por condición, 12 muestras totales, secuenciación 16S rRNA.

Este documento reúne, en orden de ejecución, todo el código usado desde la construcción de la base de datos de referencia hasta el análisis de abundancia diferencial final. Se publica junto con el manuscrito como material suplementario reproducible.

## Tabla de contenidos

1. [Entorno y dependencias](#1-entorno-y-dependencias)
2. [Construcción de la base de datos de referencia 16S](#2-construcción-de-la-base-de-datos-de-referencia-16s)
3. [Recorte de primers (cutadapt)](#3-recorte-de-primers-cutadapt)
4. [DADA2: de lecturas crudas a ASVs](#4-dada2-de-lecturas-crudas-a-asvs)
5. [Asignación taxonómica (SILVA + BLAST)](#5-asignación-taxonómica-silva--blast)
6. [Árbol filogenético (opcional)](#6-árbol-filogenético-opcional)
7. [Construcción del objeto phyloseq](#7-construcción-del-objeto-phyloseq)
8. [Filtrado de calidad y prevalencia](#8-filtrado-de-calidad-y-prevalencia)
9. [Diversidad alfa](#9-diversidad-alfa)
10. [Diversidad beta (Aitchison/CLR)](#10-diversidad-beta-aitchisonclr)
11. [Curvas de rarefacción](#11-curvas-de-rarefacción)
12. [Composición: barras, bubble plot, heatmap, Venn, UpSet](#12-composición-barras-bubble-plot-heatmap-venn-upset)
13. [Abundancia diferencial: LEfSe, ANOVA, ANCOM, ANCOM-BC2](#13-abundancia-diferencial-lefse-anova-ancom-ancom-bc2)
14. [Problemas conocidos de los paquetes y soluciones aplicadas](#14-problemas-conocidos-de-los-paquetes-y-soluciones-aplicadas)
15. [Reproducibilidad](#15-reproducibilidad)

---

## 1. Entorno y dependencias

### Herramientas de línea de comandos (conda/mamba)

```bash
mamba install -c conda-forge -c bioconda ncbi-datasets-cli entrez-direct barrnap seqkit unzip
mamba install -c bioconda cutadapt
mamba activate bioinfo
```

### Paquetes de R utilizados

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

## 2. Construcción de la base de datos de referencia 16S

Antes de poder asignar taxonomía a nivel de especie hace falta una base de referencia curada con el gen 16S de las 9 cepas concretas de SkinCom (no basta con SILVA genérico, que solo resuelve bien hasta género). `skincom_16S.sh` descarga el genoma de cada cepa desde NCBI, extrae el 16S con `barrnap`, y filtra por longitud (1300–1700 pb) para descartar fragmentos parciales.

```bash
#!/usr/bin/env bash
# =====================================================================
# skincom_16S.sh
# Obtiene el 16S rRNA de las 9 cepas de SkinCom (Lekbua et al. 2024,
# Cell Rep Methods 4:100832) a partir de sus genomas en NCBI.
#
# Dependencias (ej. con conda/mamba):
#   mamba install -c conda-forge -c bioconda ncbi-datasets-cli \
#         entrez-direct barrnap seqkit unzip
#
# Uso:
#   bash skincom_16S.sh                # todo con valores por defecto
#   THREADS=8 OUT=mi_carpeta bash skincom_16S.sh
#
# Si una cepa no se encuentra o el genoma elegido no es el correcto,
# edita OUT/cepas.tsv (4ª columna = accession manual, ej. GCF_000008345.1)
# y vuelve a lanzar el script: reutiliza lo ya descargado.
# =====================================================================
set -uo pipefail

THREADS=${THREADS:-4}
OUT=${OUT:-skincom_16S}
MINLEN=1300     # el 16S completo mide ~1500 pb; descarta hits parciales
MAXLEN=1700

for x in datasets esearch efetch esummary xtract barrnap seqkit unzip awk; do
  command -v "$x" >/dev/null 2>&1 || { echo "ERROR: falta '$x' en el PATH" >&2; exit 1; }
done

mkdir -p "$OUT"/{genomes,per_strain,fallback_nuccore,sin_match}

# ---------------------------------------------------------------------
# Tabla de cepas: etiqueta | especie | términos de cepa (separados por |) | accession manual ("-" = automático)
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
# Funciones
# ---------------------------------------------------------------------
# "A|B" -> "\"A\" OR \"B\"" para Entrez
terms_to_entrez() {
  awk -v t="$1" 'BEGIN{n=split(t,T,"|"); for(i=1;i<=n;i++) printf "%s\"%s\"", (i>1?" OR ":""), T[i]}'
}

# Lista assemblies de NCBI: acc <TAB> nivel <TAB> nombre <TAB> especie <TAB> valores infraespecíficos
list_assemblies() {  # $1 = query Entrez
  esearch -db assembly -query "$1" </dev/null 2>/dev/null \
    | esummary 2>/dev/null \
    | xtract -pattern DocumentSummary \
        -element AssemblyAccession AssemblyStatus AssemblyName SpeciesName Sub_value
}

# Elige el mejor assembly cuyo campo de cepa coincida con algún término.
# Prioridad: Complete Genome > Chromosome > Scaffold > Contig; GCF_ (RefSeq) antes que GCA_.
pick_assembly() {  # $1 = especie, $2 = términos
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

download_genome() {  # $1 = accession, $2 = carpeta destino
  local acc="$1" d="$2"
  [[ -s "$d/genome.fna" ]] && return 0
  mkdir -p "$d"
  datasets download genome accession "$acc" --include genome \
      --filename "$d/pkg.zip" --no-progressbar >/dev/null 2>&1 || return 1
  unzip -oq "$d/pkg.zip" -d "$d/pkg" || return 1
  cat "$d"/pkg/ncbi_dataset/data/"$acc"/*.fna > "$d/genome.fna" || return 1
  rm -rf "$d/pkg" "$d/pkg.zip"
}

extract_16S() {  # $1 = label, $2 = accession, $3 = carpeta genoma
  local label="$1" acc="$2" d="$3"
  barrnap --kingdom bac --threads "$THREADS" --quiet \
      --outseq "$d/rrna.fa" "$d/genome.fna" > "$d/rrna.gff" 2>/dev/null
  seqkit grep -r -n -p '^16S_rRNA' "$d/rrna.fa" 2>/dev/null \
    | seqkit seq -m "$MINLEN" -M "$MAXLEN" 2>/dev/null \
    | awk -v l="$label" -v a="$acc" \
        '/^>/{ n++; sub(/^>16S_rRNA::/,""); print ">" l "|" a "|copia" n "|" $0; next } { print }' \
    > "$OUT/per_strain/$label.16S.fa"
}

fallback_nuccore() {  # $1 = label, $2 = especie, $3 = términos
  local q; q="\"$2\"[Organism] AND ($(terms_to_entrez "$3")) AND 16S ribosomal RNA[Title] AND 1200:1700[SLEN]"
  esearch -db nuccore -query "$q" </dev/null 2>/dev/null | efetch -format fasta 2>/dev/null \
    > "$OUT/fallback_nuccore/$1.nuccore.fa"
  [[ -s "$OUT/fallback_nuccore/$1.nuccore.fa" ]]
}

# ---------------------------------------------------------------------
# Bucle principal (fd 3 para que esearch/efetch no consuman el stdin del bucle)
# ---------------------------------------------------------------------
printf 'etiqueta\tespecie\tfuente\taccession\tnivel\tn_copias_16S\tn_secuencias_unicas\tlongitud\n' \
  > "$OUT/resumen.tsv"
: > "$OUT/accesiones_usadas.tsv"

while IFS=$'\t' read -r label species terms manual <&3; do
  [[ -z "${label:-}" || "$label" == \#* ]] && continue
  echo ">>> $label ($species | $terms)"

  acc=""; level=""; fuente="genoma"
  if [[ "${manual:--}" != "-" ]]; then
    acc="$manual"; level="manual"
  else
    IFS=$'\t' read -r acc level < <(pick_assembly "$species" "$terms")
  fi

  if [[ -n "$acc" ]] && download_genome "$acc" "$OUT/genomes/$acc"; then
    echo "    genoma: $acc ($level)"
    extract_16S "$label" "$acc" "$OUT/genomes/$acc"
  else
    echo "    [!] sin genoma para la cepa; guardo candidatos de la especie en sin_match/"
    list_assemblies "\"$species\"[Organism] AND \"complete genome\"[filter] AND latest[filter]" \
      | head -n 200 > "$OUT/sin_match/$label.candidatos.tsv"
    acc="-"; level="-"; fuente="nuccore(REVISAR)"
    : > "$OUT/per_strain/$label.16S.fa"
    if fallback_nuccore "$label" "$species" "$terms"; then
      awk -v l="$label" '/^>/{ n++; sub(/^>/,""); gsub(/ /,"_"); print ">" l "|nuccore|" n "|" $0; next } { print }' \
        "$OUT/fallback_nuccore/$label.nuccore.fa" > "$OUT/per_strain/$label.16S.fa"
    else
      fuente="NINGUNA"
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
# FASTA finales
# ---------------------------------------------------------------------
cat "$OUT"/per_strain/*.16S.fa > "$OUT/SkinCom_16S_todas_las_copias.fa"
: > "$OUT/SkinCom_16S_1_por_cepa.fa"
for f in "$OUT"/per_strain/*.16S.fa; do
  [[ -s "$f" ]] && seqkit head -n 1 "$f" >> "$OUT/SkinCom_16S_1_por_cepa.fa"
done

echo
echo "=== RESUMEN ==="
column -t -s $'\t' "$OUT/resumen.tsv"
echo
echo "Salidas en $OUT/:"
echo "  SkinCom_16S_1_por_cepa.fa           (1 copia por cepa)"
echo "  SkinCom_16S_todas_las_copias.fa     (todas las copias del operón rrn)"
echo "  resumen.tsv / accesiones_usadas.tsv (para el apartado de métodos)"
echo "  sin_match/, fallback_nuccore/       (solo si alguna cepa falló: REVISAR a mano)"
```

A partir de `SkinCom_16S_1_por_cepa.fa` se construyeron dos bases derivadas: una para `dada2::assignTaxonomy()` (formato SILVA-like, hasta género) y una base BLAST (`blastdb/SkinCom_16S`) para la asignación de especie por coincidencia exacta (ver [sección 5](#5-asignación-taxonómica-silva--blast)).

---

## 3. Recorte de primers (cutadapt)

Los primers deben eliminarse antes de DADA2: su presencia interfiere con el solapamiento de lecturas, el modelo de errores y la asignación taxonómica. Región amplificada: V3–V4, primers 341F / 805R.

```bash
cd "/mnt/c/AA bioinfo tesis/comunidadmicrobianaepidermidis"

# 1. Lista de muestras
ls raw/*_1.fq.gz | xargs -n 1 basename | sed 's/_1\.fq\.gz$//' > samples.txt
cat samples.txt

# 2. Recorte de primers con cutadapt (por muestra; -g/-G = primer forward/reverse,
#    -O 15 = solapamiento mínimo, -m 200 = descarta pares que queden más cortos que
#    eso tras el recorte, --discard-untrimmed = descarta pares sin ambos primers)
cutadapt \
  -g ACTCCTACGGGAGGCAGCA \
  -G GGACTACHVGGGTWTCTAAT \
  -O 15 -m 200 --discard-untrimmed -j 0 \
  -o raw/trimmed/${sample}_R1.trimmed.fastq.gz \
  -p raw/trimmed/${sample}_R2.trimmed.fastq.gz \
  raw/${sample}_1.fq.gz raw/${sample}_2.fq.gz \
  > stats/cutadapt_primer_trimming_stats.txt
```

**Principios generales al ajustar estos parámetros para otro dataset:** comprobar la orientación de los primers, elegir los límites de longitud del amplicón apropiados para el marcador usado, e inspeccionar el informe de recorte antes de continuar — si el recorte sale mal, revisar primero la degeneración y orientación de los primers. Los parámetros a ajustar son las secuencias de primer (`-g`/`-G`, deben coincidir exactamente con el par usado) y los límites `-m`/`-M`, que deben ajustarse a la longitud esperada del amplicón tras el recorte — no copiarlos de otro proyecto sin comprobar.

**Verificación:** en este dataset, el 99.5% de los pares pasaron el filtro (152,433 → 152,433 pares con ambos primers detectados, 99.6%/99.9% R1/R2), muy por encima del umbral habitual de >90%.

---

## 4. DADA2: de lecturas crudas a ASVs

```r
library(dada2)

# Definir la ruta donde están las lecturas recortadas con Cutadapt
path <- "C:/AA Bioinfo tesis/comunidadmicrobianaepidermidis/raw/trimmed"

# Obtener rutas de los archivos R1 y R2
fnFs <- sort(list.files(path, pattern = "_R1.trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.trimmed.fastq.gz", full.names = TRUE))

# Extraer nombres de las muestras
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

# 4. Inspeccionar perfiles de calidad
## Ayudan a elegir la longitud de truncado; es uno de los pasos de criterio más
## importantes del pipeline, porque una mala elección puede destruir el
## solapamiento o dejar demasiadas bases de baja calidad.
if (!dir.exists("stats")) dir.create("stats", recursive = TRUE)

qp_F <- plotQualityProfile(fnFs[1:3])
qp_R <- plotQualityProfile(fnRs[1:3])

ggplot2::ggsave("stats/quality_profile_forward.pdf", qp_F, width = 8, height = 6)
ggplot2::ggsave("stats/quality_profile_reverse.pdf", qp_R, width = 8, height = 6)

qp_F_all <- dada2::plotQualityProfile(fnFs)
qp_R_all <- dada2::plotQualityProfile(fnRs)

ggplot2::ggsave("stats/quality_profile_forward_all.pdf", qp_F_all, width = 16, height = 12)
ggplot2::ggsave("stats/quality_profile_reverse_all.pdf", qp_R_all, width = 16, height = 12)

# 5. Filtrado y recorte
## Elimina lecturas/bases de baja calidad antes del denoising.
##
## IMPORTANTE: truncLen = c(0,0) es intencional, no un placeholder. Como cutadapt
## recorta primers tolerando indels, las lecturas resultantes tienen longitud
## VARIABLE (225-234 pb en este dataset). Un truncLen fijo (p.ej. 230) habría
## descartado cualquier lectura más corta que ese valor -- en este dataset,
## eso hubiera tirado ~78% de las lecturas. maxEE es el filtro de calidad real.
filt_path <- "C:/AA Bioinfo tesis/comunidadmicrobianaepidermidis/raw/trimmed/filtered"
dir.create(filt_path, showWarnings = FALSE)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  trimLeft = c(0, 0),
  truncLen = c(0, 0),       # CLAVE: 0 para no destruir el solapamiento tras usar Cutadapt
  maxN = 0,
  maxEE = c(2, 2),          # El filtro de calidad real: descarta lecturas con >2 errores esperados
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE,      # CLAVE: FALSE para que RStudio no se quede colgado en Windows
  matchIDs = TRUE
)

head(out)
write.csv(out, "stats/filterAndTrim_read_tracking.csv")

# 6. Aprender el modelo de error
## DADA2 modela el error de secuenciación empíricamente para distinguir
## variantes reales de errores.
set.seed(100)

errF <- learnErrors(filtFs, nbases = 1e8, multithread = FALSE, randomize = TRUE)
errR <- learnErrors(filtRs, nbases = 1e8, multithread = FALSE, randomize = TRUE)

ep_F <- plotErrors(errF, nominalQ = TRUE)
ep_R <- plotErrors(errR, nominalQ = TRUE)

ggplot2::ggsave("stats/error_rates_forward.pdf", ep_F, width = 8, height = 6)
ggplot2::ggsave("stats/error_rates_reverse.pdf", ep_R, width = 8, height = 6)

# 7. Inferencia de ASVs y ensamblaje de pares
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
head(mergers[[1]])

# 8. Tabla de secuencias y eliminación de quimeras
seqtab.all <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab.all, method = "consensus", multithread = FALSE)

dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

# Seguimiento de lecturas a través de todo el pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "sapply_seguimiento_lecturas2.csv", row.names = TRUE)
```

**Resultado:** 173 ASVs finales, retención media del 85.7% (input → no quiméricas) en las 12 muestras, sin muestras problemáticas (rango 84.8%–87.0%).

---

## 5. Asignación taxonómica (SILVA + BLAST)

Estrategia en dos pasos: `assignTaxonomy()` de DADA2 hasta género (bootstrap 50%), y BLAST contra la base curada de las 9 cepas para resolver especie por coincidencia exacta.

```r
# 1. Asignación hasta género
taxa <- assignTaxonomy(
  seqtab.nochim,
  refFasta = "refs/SkinCom_16S_dada2_taxonomy.fa",
  multithread = TRUE,
  minBoot = 50
)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

setwd("C:/AA bioinfo tesis/comunidadmicrobianaepidermidis")
library(Biostrings)

# 1) ASVs a FASTA
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV", seq_along(asv_seqs))
writeXStringSet(DNAStringSet(setNames(asv_seqs, asv_ids)), "results/asv.fa")

# 2) BLAST contra la base curada (en terminal / system2)
system2("blastn", c("-query", "results/asv.fa", "-db", "blastdb/SkinCom_16S",
                    "-outfmt", shQuote("6 qseqid sseqid pident length qcovs bitscore"),
                    "-max_target_seqs", "50", "-out", "results/asv_vs_db.tsv"))

# 3) Especie = mejor hit (los empates se reportan juntos, p.ej. "capitis/epidermidis")
hits <- read.delim("results/asv_vs_db.tsv", header = FALSE,
                   col.names = c("asv", "ref", "pident", "len", "qcov", "bitscore"))
hits$especie <- sub("_16S_[0-9]+$", "", hits$ref)

id_min <- 98; qcov_min <- 90; margen_bits <- 6     # 6 bits ~ 1 diferencia

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

# Consolidación: combinar coincidencia exacta + mejor hit BLAST en una sola
# columna Species, con el origen de la asignación registrado.
conf <- !is.na(tax_final$Species) & !is.na(tax_final$Species.1) & tax_final$Species != tax_final$Species.1
sum(conf)  # comprobación: debería ser 0 o muy pocas
tax_final[conf, c("ASV", "Genus", "Species", "Species.1")]

tax_final$Species_origen <- ifelse(!is.na(tax_final$Species), "coincidencia_exacta",
                                   ifelse(!is.na(tax_final$Species.1), "mejor_hit_BLAST", NA_character_))
tax_final$Species   <- ifelse(is.na(tax_final$Species), tax_final$Species.1, tax_final$Species)
tax_final$Species.1 <- NULL

# Géneros con una sola especie posible en la comunidad diseñada: inferir directamente
esp_unica <- c(Corynebacterium = "afermentans",
               Micrococcus     = "luteus",
               Cutibacterium   = "acnes")

falta <- is.na(tax_final$Species) & tax_final$Genus %in% names(esp_unica)
tax_final$Species[falta] <- unname(esp_unica[tax_final$Genus[falta]])
tax_final$Species_origen[falta] <- "inferida_por_genero"

table(tax_final$Species, tax_final$Species_origen, useNA = "ifany")
write.csv(tax_final, "results/taxonomy_final.csv", row.names = FALSE)
```

**Resultado:** 173 ASVs → 7 de las 9 especies del diseño recuperadas en la comunidad final (*S. hominis* y *M. luteus* ausentes — esperado, sin artefactos de mezcla entre especies). Identidad y cobertura BLAST del 100% para todas las especies asignadas por coincidencia exacta. Un pequeño número de ASVs (p.ej. ASV50, ASV72, ASV156) quedó sin resolver a nivel de especie dentro de *Staphylococcus*.

---

## 6. Árbol filogenético (opcional)

Paso opcional, usado solo si se necesita una distancia UniFrac o una figura de árbol; el análisis principal de diversidad beta de este trabajo usa distancia Aitchison/CLR (ver sección 10), que no requiere árbol.

```r
library(DECIPHER)
library(phangorn)
library(Biostrings)

seqs <- getSequences(seqtab.nochim)
asv_ids <- paste0("ASV", seq_along(seqs))     # mismos IDs que en tax_final
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

tree <- midpoint(fitGTR$tree)                 # enraizar (necesario para UniFrac)

saveRDS(tree, "results/tree.rds")
ape::write.tree(tree, "results/tree.nwk")
```

---

## 7. Construcción del objeto phyloseq

```r
# Exportar ASVs, tabla de conteos y taxonomía a los formatos que espera phyloseq
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

# Metadatos: de Excel a TSV
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

# IDs de la plataforma de secuenciación (BMK) -> nombres de condición/réplica
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

## 8. Filtrado de calidad y prevalencia

Se eliminan categorías técnicas o de muy baja información que no deben analizarse junto con la comunidad objetivo: muestras de muy baja profundidad y ASVs residuales de ruido de secuenciación.

```r
# 0) Profundidad por muestra
sort(sample_sums(ps))

# 1) Muestras con muy poca profundidad (ninguna cae con este umbral en este dataset)
ps_filtered <- prune_samples(sample_sums(ps) >= 250, ps)
nsamples(ps) - nsamples(ps_filtered)

# 2) Ruido: <0.005% de las lecturas totales o presente en una sola muestra
umbral <- 0.00005 * sum(sample_sums(ps_filtered))
otu <- as(otu_table(ps_filtered), "matrix")
if (!taxa_are_rows(ps_filtered)) otu <- t(otu)
prevalencia <- rowSums(otu > 0)
stopifnot(identical(names(prevalencia), taxa_names(ps_filtered)))

mantener <- taxa_sums(ps_filtered) >= umbral & prevalencia >= 2
ps_filtered <- prune_taxa(mantener, ps_filtered)

cat("Umbral:", round(umbral), "lecturas | ASVs:", ntaxa(ps), "->", ntaxa(ps_filtered), "\n")
cat("% lecturas conservadas:", round(100 * sum(sample_sums(ps_filtered)) / sum(sample_sums(ps)), 3), "\n")

resumen <- function(x) {
  tt <- as.data.frame(tax_table(x))
  tapply(taxa_sums(x), paste(tt$Genus, tt$Species), sum)
}
a <- resumen(ps); d <- resumen(ps_filtered)
cbind(antes = a, despues = d[names(a)])
ps_filtered
```

`ps_filtered` es el objeto usado en todos los análisis posteriores (secciones 9–13).

---

## 9. Diversidad alfa

```r
library(ggpubr)

ps_nonzero <- prune_taxa(taxa_sums(ps_filtered) > 0, ps_filtered)
sample_data(ps_nonzero)$Community <- factor(sample_data(ps_nonzero)$Community, levels = c("pHD17", "Dbap", "PF"))
table(sample_data(ps_nonzero)$Community)      # debe salir 4, 4, 4

alfa <- estimate_richness(ps_nonzero, measures = c("Observed", "Shannon", "Simpson"))
rownames(alfa) <- sample_names(ps_nonzero)       # estimate_richness cambia los guiones por puntos
alfa$Pielou <- alfa$Shannon / log(alfa$Observed)

df <- cbind(as(sample_data(ps_nonzero), "data.frame"), alfa)
df

# Con solo 4 réplicas por grupo, el test tiene poca potencia y solo detecta
# diferencias muy grandes -- ausencia de significancia no equivale a igualdad.
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
# Half-violin plots (estilo alternativo)
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

**Resultado:** sin diferencias significativas entre pHD17/Dbap/PF (Kruskal-Wallis, Shannon p=0.87) — la diversidad global de la comunidad se mantiene estable; el efecto de *bap* (ver sección 13) es específico de especie, no un cambio general en diversidad.

---

## 10. Diversidad beta (Aitchison/CLR)

Distancia Aitchison/CLR en vez de Bray-Curtis/UniFrac: más apropiada para datos composicionales de 16S (Gloor et al. 2017), y no requiere árbol filogenético.

```r
library(vegan)
library(zCompositions)

X <- as(otu_table(ps_filtered), "matrix")
if (taxa_are_rows(ps_filtered)) X <- t(X)
dim(X)
sum(X == 0)

d.czm <- cmultRepl(X, method = "CZM", label = 0)   # reemplazo de ceros
E.clr <- t(apply(d.czm, 1, function(x) log(x) - mean(log(x))))   # CLR por muestra

d.pcx    <- prcomp(E.clr)
dist.clr <- dist(E.clr)                            # distancia de Aitchison
summary(d.pcx)$importance[, 1:4]

# PCA sobre CLR (con pseudocuenta 1, dentro de phyloseq)
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

# PERMANOVA sobre las mismas distancias CLR (Aitchison)
X <- as(otu_table(ps_clr), "matrix")
if (taxa_are_rows(ps_clr)) X <- t(X)
clr_dist <- dist(X, method = "euclidean")

meta_df <- data.frame(sample_data(ps_clr))
meta_df <- meta_df[rownames(X), , drop = FALSE]
meta_df$Community <- as.factor(meta_df$Community)

set.seed(123)
adonis_res <- adonis2(clr_dist ~ Community, data = meta_df, permutations = 999)
print(adonis_res)

# Homogeneidad de dispersión entre grupos (supuesto de PERMANOVA)
anova(betadisper(clr_dist, meta_df$Community))

# ANOSIM: confirmación no paramétrica basada en rangos de distancia
set.seed(123)
ano.group <- anosim(clr_dist, meta_df$Community, permutations = 999)
ano.group

pdf("results/beta_anosim.pdf", width = 6, height = 6)
plot(ano.group)
dev.off()
```

**Resultado:** R = 0.007, P = 0.398 (ANOSIM) — las tres comunidades no se separan en el espacio de ordenación. Con n=4/grupo la potencia para separar centroides es limitada; coherente con el resto del análisis (el efecto de *bap* es específico de una especie, no un desplazamiento composicional global).

---

## 11. Curvas de rarefacción

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
  theme_classic() + xlab("Lecturas") + ylab("Shannon")

print(p)
ggsave("stats/shannon_rarefaction.pdf", p, width = 7, height = 5)
```

---

## 12. Composición: barras, bubble plot, heatmap, Venn, UpSet

### 12.1 Barras de abundancia relativa (por muestra y por comunidad)

```r
library(dplyr)
library(ggplot2)

nivel <- "Species"

ps_rel <- ps_filtered %>%
  tax_glom(taxrank = nivel, NArm = FALSE) %>%
  transform_sample_counts(function(x) x / sum(x))

# Por muestra
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

# Por comunidad (media de las réplicas)
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

# Tabla de porcentajes (las barras no dejan ver bien las especies minoritarias)
tt <- as.data.frame(type_tax); tt$taxon <- tt[[nivel]]
round(100 * xtabs(Abundance ~ taxon + Community, data = tt), 2)
```

### 12.2 Bubble plot (media ± SD por especie y comunidad)

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

# --- 1) SD entre réplicas, por especie y comunidad --------------
# Calculado ANTES de fusionar muestras: cada réplica se normaliza a
# % por su cuenta, y luego se calcula la SD dentro de cada Community.
sd1 <- ps_filtered %>%
  tax_glom(taxrank = nivel, NArm = FALSE) %>%
  transform_sample_counts(function(x) x / sum(x) * 100) %>%
  psmelt() %>%
  group_by(OTU, Community) %>%
  summarize(SD = sd(Abundance, na.rm = TRUE), .groups = "drop")

# --- 2) Media por comunidad --------------------------------------
bub1 <- merge_samples(ps_filtered, "Community")

type_taxa1 <- bub1 %>%
  tax_glom(taxrank = nivel, NArm = FALSE) %>%
  transform_sample_counts(function(x) x / sum(x) * 100) %>%
  psmelt() %>%
  mutate(Community = as.character(Sample)) %>%   # merge_samples desordena metadatos; reconstruir Community
  arrange(OTU)

dat1 <- type_taxa1 %>%
  left_join(sd1, by = c("OTU", "Community")) %>%
  as.data.table()

# --- 3) Especies raras NO se colapsan -----------------------------
# A diferencia de una comunidad ambiental ruidosa con cientos de taxones,
# esta es una comunidad sintética DEFINIDA de 9 especies: si una especie
# cae mucho en una condición (p.ej. el fondo Dbap), eso es justo lo que
# se quiere ver, no ocultar. Solo se colapsa el ~0.01% de lecturas sin
# asignación de especie.
dat1[, Species_plot := as.character(get(nivel))]
dat1[is.na(Species_plot) | Species_plot == "NA", Species_plot := "Unassigned (background noise, ~0.01%)"]

# Orden del eje Y Y de la leyenda: por abundancia media descendente.
# Asignar los MISMOS niveles de factor antes de graficar mantiene el
# orden de color de la leyenda sincronizado con las burbujas del eje Y.
orden_especies <- dat1[, .(m = mean(Abundance, na.rm = TRUE)), by = Species_plot][order(-m), Species_plot]
dat1[, Species_plot := factor(Species_plot, levels = orden_especies)]
sample_order <- c("pHD17", "Dbap", "PF")

# --- 4) Gráfico -----------------------------------------------------
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
  tax_aggregate = "Species",   # agrega a especie, no OTU
  tax_add = "Genus",           # añade el género junto al nombre de especie en el eje Y
  tax_show = n_species,        # muestra todas las especies definidas, ninguna oculta
  color_vector = c("snow", "steelblue"),
  plot_values = TRUE
) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Most abundant species per community")

print(xx)
ggplot2::ggsave("results/heatmap_especies.pdf", xx, width = 5, height = 6)
```

### 12.4 Diagrama de Venn (ampvis2)

```r
# ================================================================
# Venn diagram of species shared/unique across communities
# ================================================================
library(ampvis2)
library(phyloseq)

# amp_venn no tiene parámetro de agregación como amp_heatmap's tax_aggregate,
# así que se agrega a especie primero (tax_glom), igual que en bubble/heatmap.
ps_species <- tax_glom(ps_filtered, taxrank = "Species", NArm = FALSE)
ampvis2_species <- amp_load(ps_species)

venn <- amp_venn(
  ampvis2_species,
  group = "Community",
  cut_a = 0.01,   # abundancia relativa media mínima (%) para contar una especie como "presente"
  cut_f = 75,     # % mínimo de réplicas del grupo donde debe aparecer (75% = 3 de 4)
  detailed_output = TRUE
)

venn$plot
```

**Resultado:** las 6 especies con abundancia >0.01% son compartidas por las 3 condiciones (core del 100%) — ninguna especie exclusiva de una sola condición.

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

Con solo 3 grupos y un core del 100% (ver Venn), el UpSet plot no aporta información adicional respecto al Venn diagram — se incluye aquí por completitud pero se descartó de la presentación final de resultados.

---

## 13. Abundancia diferencial: LEfSe, ANOVA, ANCOM, ANCOM-BC2

Cuatro métodos independientes, con lógicas estadísticas distintas, para evaluar si el gen *bap* funcional confiere una ventaja a alguna especie de la comunidad. Los cuatro se ejecutan sobre `ps_filtered` agregado a nivel de especie (`tax_glom`).

### 13.1 LEfSe

> **Nota de reproducibilidad:** `run_lefse()` de `microbiomeMarker` (v1.13.2, Bioconductor devel) tiene un bug real confirmado en su propio código fuente — ver [sección 14](#14-problemas-conocidos-de-los-paquetes-y-soluciones-aplicadas). El código de abajo usa una copia parcheada, `run_lefse_fixed()`.

```r
# ================================================================
# LEfSe across the three communities, with microbiomeMarker
# ================================================================
library(microbiomeMarker)
library(phyloseq)

# PATCH: ver sección 14 para la explicación completa del bug y el fix.
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
    sig_otus <- otus_test[, sig_ind, drop = FALSE]   # <-- EL FIX (era: otus_test[, sig_ind])
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

**Resultado:** 1 marcador — *S. epidermidis*, enriquecida en pHD17 (Kruskal-Wallis p≤0.05, LDA≥2).

### 13.2 ANOVA (paramétrico sobre CLR)

> **Nota de reproducibilidad:** `run_test_multiple_groups()` de `microbiomeMarker` se descartó aquí — su parámetro `transform` no acepta `"clr"`, y enruta a través del mismo `preprocess_ps()` interno que rompe con datos ya CLR-transformados (ver sección 14). En su lugar: ANOVA de una vía por especie, directo y transparente con `stats::aov()`.

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

**Resultado:** ninguna especie sobrevive corrección FDR (q<0.05); *S. epidermidis* es la más cercana (F=5.85, p=0.024, q=0.165).

### 13.3 ANCOM

```r
# ================================================================
# ANCOM (original, Mandal et al. 2015)
# run_ancom() hace su PROPIO manejo composicional interno (confirmado
# vía args(run_ancom): transform solo acepta identity/log10/log10p,
# sin "clr", y tiene su propio norm="TSS" por defecto) -- se le pasa
# ps_species crudo, no el objeto CLR-transformado.
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
ancom_df <- as(ancom, "data.frame")   # necesario: write.csv() directo sobre el objeto S4 falla (ver sección 14)
write.csv(ancom_df, "results/markers_ANCOM.csv", row.names = FALSE)
```

**Resultado:** ningún marcador identificado al umbral estándar (W_cutoff=0.75). Con solo 7 especies, ANCOM exige que una especie sea significativa frente al 75% de las demás en comparación por pares — un listón muy alto para tan pocos taxones. Coherente con el ANOVA-FDR.

### 13.4 ANCOM-BC2 (modelo combinado, 3 grupos)

```r
# ================================================================
# ANCOM-BC2 across the three communities
# Sin covariable real disponible en metadata.tsv (solo Community, un
# id, y DNA ng/uL -- valor técnico de extracción, no incluido).
# Referencia = Dbap: mismo plásmido que pHD17, bap mutado -- aísla el
# efecto de bap (pHD17 vs Dbap) del efecto de carga del plásmido
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
  prv_cut = 0.01,       # laxo a propósito: solo 7-9 especies diseñadas, no filtrar ninguna por accidente
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
  trend = FALSE   # sin orden natural entre pHD17/Dbap/PF; trend=TRUE exige matrices de contraste que no hacen falta
)

res <- output_combined$res
res
write.csv(res, "results/markers_ANCOMBC2.csv", row.names = FALSE)
```

> **Nota de reproducibilidad importante:** con `pairwise=TRUE` + `dunnet=TRUE` a la vez, `output_combined$res_pair` tiene las columnas `p_`/`q_` rotas (ver sección 14) — para las comparaciones por pares, usar la versión de la sección 13.5 en su lugar. `res` (sin pairwise) sí es de fiar para los dos contrastes de referencia (pHD17 vs Dbap, PF vs Dbap).

### 13.5 ANCOM-BC2 (3 comparaciones por pares, limpias)

```r
# ================================================================
# ANCOM-BC2, como 3 comparaciones de dos grupos independientes en vez
# de la ejecución combinada pairwise+dunnet de la sección 13.4.
#
# Por qué: con global/pairwise/dunnet=TRUE a la vez, las columnas p/q
# de res_pair salen rotas (p=1 para todo, incluso un caso con W=3.4
# que debería dar p=0.009 -- bug real confirmado, no nuestro). Sub-
# setear a dos grupos y correr ancombc2 simple (pairwise=FALSE etc.)
# es la solución estándar para factores multinivel y evita por
# completo ese camino de código roto.
# ================================================================
library(ANCOMBC)
library(phyloseq)
library(dplyr)

run_pair <- function(ps, level_ref, level_other, contrast_label) {
  # prune_samples(), no subset_samples(): subset_samples() usa NSE
  # internamente y no encuentra level_ref/level_other cuando se llama
  # dentro de una función (limitación real de phyloseq, no nuestro bug)
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

**Resultado:** pHD17 vs Dbap: lfc=1.24, q=0.095 (la más cercana a significancia, *S. epidermidis*). PF vs Dbap y pHD17 vs PF: sin marcadores significativos ni sugestivos.

### 13.6 Visualización de ANCOM-BC2 (bar plot + volcano, 3 contrastes)

```r
# ================================================================
# ANCOM-BC2 visualizations: all 3 pairwise contrasts
# Requiere haber corrido la sección 13.5 antes -- usa su `df_res`.
# ================================================================
library(dplyr)
library(ggplot2)

orden_taxa <- df_res %>%
  dplyr::filter(contrast == "pHD17 vs Dbap") %>%
  dplyr::arrange(lfc) %>%
  dplyr::pull(taxon)
df_res$taxon <- factor(df_res$taxon, levels = orden_taxa)

# ---- 1. Bar plot de log fold change, 3 contrastes ----
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

# ---- 2. Volcano plot, 3 contrastes ----
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

### 13.7 Síntesis de los 4 métodos

| Método | ¿Corrección múltiple? | Resultado para *S. epidermidis* |
|---|---|---|
| LEfSe | No (KW + Wilcoxon + LDA) | Único marcador, enriquecida en pHD17 |
| ANOVA (CLR) | Sí (FDR) | p=0.024, q=0.165 — más cercana, no significativa |
| ANCOM | Implícita (W_cutoff) | Ningún marcador al umbral estándar |
| ANCOM-BC2 (pares) | Sí (FDR) | pHD17 vs Dbap: q=0.095 — la más cercana |

Validado además de forma independiente contra el informe comercial de secuenciación (Metastats/BMK): p crudo <0.05 en las 3 comparaciones por pares para *S. epidermidis*, pero q>0.05 en las 3 — mismo patrón exacto, confirmando que ninguno de los dos análisis se contradice, solo difieren en si destacan el p crudo o el q corregido como criterio de "significativo".

**Conclusión:** ningún método alcanza significancia formal tras corrección por comparaciones múltiples (n=4/grupo limita la potencia), pero la convergencia direccional entre métodos con lógicas estadísticas distintas señala consistentemente a *S. epidermidis* como la especie más afectada por la funcionalidad de *bap*.

---

## 14. Problemas conocidos de los paquetes y soluciones aplicadas

Documentado aquí porque son bugs reales de las versiones instaladas (no errores de uso), reproducibles y relevantes para cualquiera que reutilice este código con las mismas versiones de paquete.

### 14.1 `microbiomeMarker::run_lefse()` — vector desnudo al seleccionar una sola columna

**Paquete/versión:** `microbiomeMarker` 1.13.2 (Bioconductor devel).

**Síntoma:** `Error: Can't recycle .x (size N) to match .y (size 0)` dentro de `purrr::map2_lgl()`, de forma intermitente según cuántos taxones pasan el filtro de Kruskal-Wallis.

**Causa raíz:** en el código fuente de `run_lefse()`, la línea `sig_otus <- otus_test[, sig_ind]` usa el subsetting por defecto de `data.frame` de R (`drop = TRUE`), que convierte el resultado a un vector plano cuando exactamente **una** columna pasa el filtro — rompiendo la siguiente llamada a `purrr::map2_lgl()`, que espera una lista/data.frame.

**Solución:** copia exacta de la función con `drop = FALSE` añadido en esa línea (`run_lefse_fixed()`, sección 13.1), con `environment()` reasignado al namespace del paquete para que siga viendo sus funciones internas no exportadas (`preprocess_ps`, `test_rep_wilcoxon`, etc.).

### 14.2 `preprocess_ps()` interno — malinterpreta datos ya CLR-transformados

**Síntoma:** `The library size of sample(s): X is/are zero` al pasar un objeto ya transformado a CLR (valores negativos) a `run_lefse()` o `run_test_multiple_groups()`.

**Causa raíz:** el helper interno `preprocess_ps()`, compartido por varias funciones `run_*` del paquete, calcula el "tamaño de librería" como la suma de valores por muestra — una operación que asume conteos crudos no negativos. Con datos CLR (media ~0 por muestra, mitad de los valores negativos), esa suma da un número cercano a cero, y la función lo interpreta como profundidad de secuenciación nula.

**Solución:** nunca pre-transformar a CLR antes de pasar los datos a estas funciones — dejar que su propio parámetro `norm`/`transform` haga la normalización internamente (`run_lefse(norm="CPM")`, `run_ancom(norm="TSS")`). Para el ANOVA, evitar la función del paquete por completo y usar `stats::aov()` directo sobre la matriz CLR (sección 13.2).

### 14.3 `ancombc2()` — `p`/`q` degenerados en `res_pair` con `pairwise` + `dunnet` combinados

**Paquete/versión:** `ANCOMBC` (release de Bioconductor, no devel — a diferencia de los dos anteriores).

**Síntoma:** con `pairwise=TRUE` y `dunnet=TRUE` simultáneamente, todas las columnas `p_`/`q_` de `output_combined$res_pair` salen exactamente en `1`, sin importar el valor real del estadístico `W`.

**Causa raíz confirmada:** comparando fila a fila `res` (modelo simple, sin pairwise) contra `res_pair` para el mismo coeficiente nominal, `lfc`, `se` y `W` son **idénticos** entre ambos objetos (p.ej. W=3.397898 en los dos), pero `p` pasa de 0.0094 (correcto) a 1 (imposible para ese W). El modelo subyacente es el mismo; algo en el cálculo conjunto de p-valores de la rama `pairwise+dunnet` con n pequeño (n=4/grupo aquí) devuelve un valor degenerado.

**Solución:** no usar `res_pair` de una ejecución combinada. En su lugar, para cada par de grupos, ejecutar `ancombc2()` de forma independiente y simple (`pairwise=FALSE`, `dunnet=FALSE`, `global=FALSE`) sobre un subset de solo esos dos grupos (`prune_samples()`, sección 13.5) — exactamente tan fiable como la ejecución de referencia que sí funcionó bien.

### 14.4 `phyloseq::subset_samples()` — evaluación no estándar rota dentro de una función

**Síntoma:** `error al evaluar el argumento 'table' ... objeto 'level_ref' no encontrado`, al llamar `subset_samples(ps, Community %in% c(level_ref, level_other))` dentro de una función definida por el usuario.

**Causa raíz:** `subset_samples()` usa evaluación no estándar (NSE) internamente y no resuelve correctamente variables definidas en el entorno de la función que la llama — una limitación conocida de phyloseq, no específica de este análisis.

**Solución:** usar `prune_samples()` con un vector lógico calculado explícitamente (`sample_data(ps)$Community %in% c(...)`) en vez de `subset_samples()`, siempre que el subsetting ocurra dentro de una función.

### 14.5 `write.csv()` sobre un objeto S4 `marker_table` con filas

**Síntoma:** `Error en validObject(.Object): invalid class "marker_table" object: marker table must have non-zero dimensions`, al llamar `write.csv()` directamente sobre el resultado de `marker_table()` cuando SÍ hay marcadores (con cero marcadores, `marker_table()` devuelve `NULL` y no hay error, distinto comportamiento).

**Causa raíz:** `marker_table` es una clase S4; `write.csv()` intenta coaccionarla internamente a `data.frame`, y esa ruta de coerción vuelve a disparar la validación de la clase S4 `marker_table`, que falla en ese punto del pipeline de `write.table()`.

**Solución:** convertir explícitamente con `as(x, "data.frame")` antes de `write.csv()` (sección 13.3) — evita por completo la ruta de coerción interna problemática.

---

## 15. Reproducibilidad

Pendiente de incluir en la versión final publicada junto con el manuscrito:

```r
sessionInfo()
# o, si el proyecto usa renv:
# renv::snapshot()
```

**Resumen de versiones relevantes usadas en este análisis:**
- `microbiomeMarker` 1.13.2 (Bioconductor devel) — ver bugs documentados en 14.1–14.2
- `ANCOMBC` (release Bioconductor) — ver bug documentado en 14.3
- `phyloseq`, `dada2`, `vegan`, `ampvis2`, `UpSetR`, `zCompositions`, `compositions`, `microbiome` — versión exacta pendiente de `sessionInfo()`

**Datos de entrada:** 12 muestras (4 réplicas × 3 condiciones: pHD17, Dbap, PF), secuenciación 16S V3–V4 (BMK, proyecto BMKO260810-AA254-01-01). Base de datos de referencia: 9 cepas de SkinCom (Lekbua et al. 2024).
