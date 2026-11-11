# Module 3: R for prevalence analysis

R is the standard language for statistics and figures in biology. In this module you go from a presence/absence table to prevalences with confidence intervals, hypothesis tests and publication-ready plots.

[Back to the course map](../README.md)

- [Lesson 14: R basics](#lesson-14-r-basics)
- [Lesson 15: Prevalence with the tidyverse](#lesson-15-prevalence-with-the-tidyverse)
- [Lesson 16: Statistics for prevalence data](#lesson-16-statistics-for-prevalence-data)
- [Lesson 17: Plots with ggplot2](#lesson-17-plots-with-ggplot2)

---

## Setup

1. Install **R** from https://cran.r-project.org. On Apple Silicon choose the `arm64` package.
2. Install **RStudio Desktop** from https://posit.co/download/rstudio-desktop/.
3. In RStudio: **File > New Project > Existing Directory** and choose `~/learn`. Working inside a project means paths such as `data/strains.csv` work without changing directories.
4. Install the packages you need (once):

```r
install.packages("tidyverse")
```

In RStudio, write code in a script (**File > New File > R Script**) and run the current line with **Cmd + Enter**. Save scripts in `scripts/` with the `.R` extension.

---

## Lesson 14: R basics

**Goal:** understand values, vectors and data frames.

### Values and variables

```r
2 + 3
n_strains <- 10
n_positive <- 6
100 * n_positive / n_strains
```

`<-` assigns a value to a name (read it as "gets"). Comments start with `#`.

### Vectors

A vector is a list of values of the same type:

```r
bap <- c(1, 1, 0, 0, 1, 1, 0, 1, 0, 1)
length(bap)
sum(bap)
mean(bap) * 100
```

Because presence is coded as 0/1, the **mean is the prevalence**.

Indexing and logical tests:

```r
bap[1]
bap[1:3]
bap == 1
which(bap == 1)
```

### Data frames

A data frame is a table: columns are vectors of the same length.

```r
strains <- read.csv("data/strains.csv")
head(strains)
str(strains)
nrow(strains)
```

Columns are accessed with `$`:

```r
strains$bap
mean(strains$bap) * 100
table(strains$origin)
table(strains$origin, strains$bap)
```

Filter rows with `[rows, columns]`:

```r
strains[strains$origin == "clinical", ]
strains[strains$bap == 1, c("strain", "origin")]
```

### Getting help

```r
?mean
```

### Exercise 14

1. Create a vector with the *esp* values of the 10 strains and compute its prevalence.
2. Using `strains`, compute the prevalence of *bap* in clinical strains only.
3. What does `table(strains$bap, strains$esp)` show?

---

## Lesson 15: Prevalence with the tidyverse

**Goal:** compute prevalences by group in a readable way.

```r
library(tidyverse)
strains <- read_csv("data/strains.csv")
```

### The pipe `|>`

`|>` passes the result on the left into the function on the right, like `|` in the Terminal. Read it as "then".

### Core verbs

| Verb | What it does |
|---|---|
| `filter()` | Keeps rows that meet a condition |
| `select()` | Keeps columns |
| `mutate()` | Creates or changes columns |
| `group_by()` | Defines groups |
| `summarise()` | Reduces each group to one row |
| `arrange()` | Sorts rows |

```r
strains |> filter(origin == "clinical")

strains |>
  summarise(n = n(), positives = sum(bap), prevalence = 100 * mean(bap))

strains |>
  group_by(origin) |>
  summarise(n = n(),
            bap_prev = 100 * mean(bap),
            esp_prev = 100 * mean(esp))
```

### Long format: one row per strain and gene

With many genes, it is easier to have one row per strain-gene pair:

```r
long <- strains |>
  pivot_longer(c(bap, esp), names_to = "gene", values_to = "present")

long
```

Now one summary works for every gene:

```r
prev <- long |>
  group_by(origin, gene) |>
  summarise(n = n(),
            positives = sum(present),
            prevalence = 100 * positives / n,
            .groups = "drop")

prev
write_csv(prev, "results/prevalence_by_origin.csv")
```

With 30 genes, you would select them with `pivot_longer(-c(strain, origin), ...)` instead of naming each one.

### Joining tables

The presence matrix from Lesson 13 and your metadata share the `strain` column:

```r
matrix <- read_csv("results/presence_matrix.csv")
meta <- read_csv("metadata.csv")

data <- meta |> left_join(matrix, by = "strain")
```

`left_join` keeps every strain in `meta`. Check for strains with missing gene data:

```r
data |> filter(is.na(NUM_FOUND))
```

### Exercise 15

1. Add a column `both` that is 1 when a strain carries both genes: `mutate(both = as.integer(bap == 1 & esp == 1))`.
2. Compute the prevalence of `both` by origin.

---

## Lesson 16: Statistics for prevalence data

**Goal:** report prevalences with confidence intervals and test whether they differ between groups.

### Confidence intervals

A prevalence of 5/6 (83%) is very uncertain with only 6 strains. Always report the interval:

```r
binom.test(5, 6)$conf.int
```

This is the exact Clopper-Pearson interval. The Wilson interval is a common alternative:

```r
prop.test(5, 6, correct = FALSE)$conf.int
```

For every group and gene at once:

```r
prev_ci <- prev |>
  rowwise() |>
  mutate(ci_low  = 100 * binom.test(positives, n)$conf.int[1],
         ci_high = 100 * binom.test(positives, n)$conf.int[2]) |>
  ungroup()

prev_ci
```

`rowwise()` makes `mutate` work one row at a time, which `binom.test` needs.

### Comparing two groups: Fisher's exact test

Is *bap* more frequent in clinical than in commensal strains?

```r
tab <- table(strains$origin, strains$bap)
tab
fisher.test(tab)
```

- The **p-value** is the probability of a difference at least this large if there were no real association.
- The **odds ratio** measures the strength of the association. Its direction depends on the order of rows and columns in the table, so check `tab` before interpreting it.
- Use Fisher's test rather than the chi-squared test when any expected count is below 5, which is usual with clinical collections.

### Many genes: correct for multiple testing

Testing 20 genes at p < 0.05 gives about one false positive by chance. Adjust the p-values:

```r
tests <- long |>
  group_by(gene) |>
  filter(sum(present) > 0, sum(present) < n()) |>
  summarise(p_value = fisher.test(table(origin, present))$p.value) |>
  mutate(p_adj = p.adjust(p_value, method = "BH"))

tests
```

- The `filter` drops genes present in all or no strains, which cannot be tested.
- `BH` (Benjamini-Hochberg) controls the false discovery rate and is the usual choice.

### Interpretation checklist

- Report n, positives, prevalence and 95% CI for every group.
- Strains from the same outbreak or clone are not independent; many related isolates can create false associations. Describe your sampling and, if possible, account for lineage (for example MLST type).
- Association is not causation: a gene enriched in clinical strains is a hypothesis to test in the lab.

### Exercise 16

1. Compute the 95% CI of *bap* prevalence in clinical strains (5/6). How wide is it?
2. Run Fisher's test for *bap* by origin. Is the difference significant? Why might that be, even though 83% vs 25% looks large?

---

## Lesson 17: Plots with ggplot2

**Goal:** make clear figures for the thesis.

A ggplot is built in layers: data, aesthetics (`aes`: what goes on each axis or colour), geometries (`geom_`: bars, points, tiles) and formatting.

### Prevalence with confidence intervals

```r
ggplot(prev_ci, aes(x = gene, y = prevalence, fill = origin)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Gene", y = "Prevalence (%)", fill = "Origin") +
  theme_classic()

ggsave("results/prevalence.png", width = 6, height = 4, dpi = 300)
```

`ggsave` saves the last plot. Use `.pdf` instead of `.png` for a vector figure.

### Presence/absence heatmap

```r
ggplot(long, aes(x = gene, y = strain, fill = factor(present))) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c("0" = "grey90", "1" = "steelblue"),
                    labels = c("Absent", "Present"), name = NULL) +
  facet_grid(origin ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_minimal()

ggsave("results/heatmap.png", width = 4, height = 6, dpi = 300)
```

`facet_grid(origin ~ .)` splits the strains into one panel per origin.

### Tips for thesis figures

- Label axes with units, for example "Prevalence (%)".
- Use the same colour for the same group in every figure.
- Prefer `theme_classic()` or `theme_minimal()` over the default grey background.
- Show n in the figure or the legend.

### Exercise 17

1. Plot only *esp*, with origin on the x axis.
2. Save it as a PDF of 5 x 4 inches.

Next: [Module 4, Python](module-4-python.md)
