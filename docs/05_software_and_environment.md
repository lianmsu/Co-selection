# Software and environment

Versions are the ones reported in the manuscript (Methods 2.3–2.4); the
production runs were performed with the conda environments named below.

## Bioinformatics software

| Step | Software | Version | Where |
| --- | --- | --- | --- |
| Amplicon merging, denoising, ASV table, taxonomy | usearch | v10.0.240 | `pipelines/01` |
| Read QC, assembly, binning, bin refinement | metaWRAP | v1.2.1 | not included (standard commands) |
| Binning (inside metaWRAP) | CONCOCT / MetaBAT2 / MaxBin2 | v1.0.0 / v2.12.1 / v2.2.6 | not included |
| MAG dereplication | dRep | v3.4.2 | not included |
| MAG taxonomy and tree | GTDB-Tk | v2.3 (GTDB R214) | `pipelines/02` |
| ORF prediction | Prodigal | v2.6.3 | not included |
| Protein search against SARG/BacMet/AromaDeg/MGE | diamond | v2.0.6 | `pipelines/03` |
| Read-based gene quantification | ARGs-OAP | v3.2.2 | `pipelines/04` |
| Plasmid detection on MAGs | PlasFlow | v1.1 | not included |
| Network visualisation | Gephi | v0.10.1 | after `analyses/03` |
| Tree annotation | iTOL | web | `analyses/02`, section 1 |

## Conda environments

```
gtdbtk2.3     GTDB-Tk 2.3 + pplacer + mash database (GTDB R214)
args_oap      ARGs-OAP v3.2.2 (stage_one / stage_two)
```

The remaining steps used the software as installed in the working
environment; only the commands are given.

## R

The manuscript reports R v4.3.1. The scripts were checked to parse under
R 4.6.0 and use:

```
ape            corrplot       data.table     dplyr          EcoSimR
ggbreak        ggcorrplot     ggpmisc        ggpubr         ggsignif
ggh4x          ggplot2        ggvenn         gtable         Hmisc
igraph         itol.toolkit   MASS           minpack.lm     patchwork
pheatmap       randomForest   RColorBrewer   readxl         reshape2
rfPermute      rstatix        scales         splines        stringr
tidyr          tidyverse      vegan          wesanderson
```

Install with, for example:

```r
install.packages(c("ape", "corrplot", "data.table", "dplyr", "ggpubr",
                   "ggsignif", "ggh4x", "ggplot2", "ggvenn", "ggbreak",
                   "ggcorrplot", "ggpmisc", "gtable", "Hmisc", "igraph",
                   "MASS", "minpack.lm", "patchwork", "pheatmap",
                   "randomForest", "RColorBrewer", "readxl", "reshape2",
                   "rstatix", "scales", "stringr", "tidyr", "tidyverse",
                   "vegan", "wesanderson"))
# BiocManager-style / GitHub packages
# EcoSimR, rfPermute, itol.toolkit (see the package documentation)
```

## Databases

| Database | Used for | Version / file |
| --- | --- | --- |
| SARG | ARGs (read-based quantification and diamond annotation) | v3.2 (`4.SARG_v3.2_20220917_Short_subdatabase`, `SARG_v3.2_diamond`) |
| BacMet | MRGs | as distributed with the diamond database |
| AromaDeg | ADGs (aromatic-compound degradation) | `AromaDeg.fasta` |
| MGE database | MGEs (integrase / transposase / plasmid markers) | `MGE.fasta` |
| MiDAS | 16S taxonomy (SINTAX) | 4.8.1 |
| GTDB | MAG taxonomy (GTDB-Tk) | R214 |

## Notes and small fixes

* **Primer trimming.** The manuscript text quotes `-stripleft 25
  -stripright 20`; the script executed for the study uses 27 and 28. The
  script values are kept as the defaults and can be overridden with
  `PREFIX_TRIM` / `SUFFIX_TRIM`.
* **usearch binary.** The original script called `usearch11` in the SINTAX
  step although the rest of the workflow uses usearch 10; the wrapper now uses
  a single `USEARCH` variable.
* **`turn_to_table.py` shebang.** It pointed to a hard-coded conda interpreter;
  it now uses `/usr/bin/env python3` and `os.path.basename`, which also makes it
  work with Windows paths.
* **`gene_abundance.sh`.** The four databases were run by four copies of the
  same two commands; they are now one `run_args_oap` function called four
  times, with identical parameters.
* **`diamond_res.sh`.** Split into `04_diamond_results_to_tables.sh`
  (merging + structure annotation) and `05_gene_counts_per_genome.sh`
  (per-genome counts). The `tongji2.txt` header is now written with `>`
  instead of `>>` so that reruns do not append a second header.
