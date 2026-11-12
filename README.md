# Lab-22

Code, pipelines and training material from my PhD thesis at IdAB-CSIC on gut pathobiont bacteria that produce amyloid-forming surface proteins (Bap, Esp).

## Contents

| Folder / file | What it is |
|---|---|
| [`bioinfo-course/`](bioinfo-course/README.md) | **Bioinformatics from Scratch**: a hands-on course for life scientists who have never programmed. It goes from opening the Terminal to a reproducible analysis of `bap` and `esp` prevalence in clinical genomes (bash, conda, BLAST, abricate, R, Python, Git). |
| [`Microbiome_Analysis.md`](Microbiome_Analysis.md) | 16S rRNA analysis of a 9-species synthetic community (pHD17, Dbap and PF conditions): reference database, cutadapt, DADA2 and differential abundance. |
| [`rnaseq_F1/`](rnaseq_F1/) | Python RNA-seq pipeline: SRA download → bowtie2 alignment → sorted BAM → strand-specific WIG tracks, for paired-end and single-end reads. |
| [`gel_analysis/`](gel_analysis/) | `measure_gel.py`: detects wells and bands in an agarose gel photo and estimates band sizes from the ladder (log-linear fit, reference bands anchored by brightness). Can also draw an annotated figure (lane numbers, size scale, arrows) ready for a presentation. Needs numpy and Pillow. |

## Repository structure

```
.
├── README.md
├── Microbiome_Analysis.md      # 16S synthetic community pipeline
├── rnaseq_F1/                  # RNA-seq pipeline (Python)
│   ├── pipeline_paired_end.py
│   ├── pipeline_single_end.py
│   └── wig_utils.py            # helper functions imported by both pipelines
├── gel_analysis/
│   └── measure_gel.py          # agarose gel band sizing from a photo
└── bioinfo-course/             # introductory bioinformatics course
    ├── README.md               # course map: start here
    ├── lessons/                # modules 0-5 and the final project
    ├── data/                   # small made-up example files
    ├── cheatsheet.md
    └── solutions.md
```
