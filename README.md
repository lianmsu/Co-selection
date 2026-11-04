# Co-selection

Code accompanying **"Understanding pollutant-driven shifts of antibiotic resistome in
activated sludge: A lab-scale study"** — Li P., Yin Z., Ye L.
*Journal of Hazardous Materials* **487** (2025) 137199.
[doi:10.1016/j.jhazmat.2025.137199](https://doi.org/10.1016/j.jhazmat.2025.137199)

The repository contains every script used to go from raw sequencing reads to the figures of
the manuscript: the usearch pipeline for the 16S rRNA gene amplicon data, the GTDB-Tk, diamond
and ARGs-OAP pipelines for the metagenomes, and the three R notebooks that produce the
statistics and the figures.

> 中文版说明见 **[README_zh.md](README_zh.md)**。

---

## Study in brief

Twenty-one lab-scale bioreactors were operated with activated sludge from an industrial
wastewater treatment plant and exposed for 30 days to three heavy metals (Cr(VI), Cd(II),
Cu(II)) or three aromatic organics (phenol, 4-aminophenol, 4-nitrophenol); a further 30 days
followed without pollutant dosing. Fifteen samples (14 sludge samples plus the seed sludge)
were analysed by 16S rRNA gene amplicon and shotgun metagenomic sequencing.

* Both pollutant classes increased the total ARG abundance (1.7–2.9-fold over the seed
  sludge), with `bacA` and `sul1` the most abundant ARGs throughout.
* ARG composition separated heavy-metal from organic treatments (NMDS, stress = 0.11),
  and the elevated ARG abundance tracked the pollutant-driven shift of the bacterial
  community rather than the pollutants themselves.
* Co-occurrence and Mantel analyses point to an **indirect co-selection** of ARGs through
  metal resistance genes (MRGs) and aromatic-compound degradation genes (ADGs), concentrated
  in Gammaproteobacteria.
* Although ARG abundance rose, the proportion of high-risk ARGs did not, i.e. abundance alone
  does not equal health risk.

## What is in here

```
Co-selection/
├── pipelines/          # upstream bioinformatics (shell + one python helper)
│   ├── 01_16S_amplicon/                usearch: merge, trim, UNOISE3, ASV table, SINTAX
│   ├── 02_metagenome_binning/          GTDB-Tk taxonomy and tree of the 257 MAGs (Fig. 3d)
│   ├── 03_resistance_gene_annotation/  diamond vs SARG/MRG/ADG/MGE, tables per MAG (Fig. 3d, 3e)
│   └── 04_gene_abundance/              ARGs-OAP v3.2.2 quantification (Fig. 1, 2, 5)
├── analyses/           # downstream statistics and figures (R)
│   ├── 01_16S_community/               Fig. 6, Figs. S7, S8
│   ├── 02_metagenome_resistome/        Fig. 1, 2, 3c, 5, Figs. S2-S6
│   └── 03_MAG_gene_annotation/         Fig. 3d, Fig. 4, Fig. S4
├── data/               # small input tables (see data/README.md)
└── docs/               # workflow, layout, file mapping, inputs/outputs, software
```

Both layers are numbered in execution order; each file carries a header with its purpose, the
manuscript figure it produces, its main outputs and its original path in the working copy.
The two long R scripts were single sequential notebooks in the original repository, so they
also carry a section index; [`docs/04`](docs/04_script_inputs_outputs.md) expands it into a table.

## Documentation

| Document | Content |
| --- | --- |
| [docs/01_pipeline_overview.md](docs/01_pipeline_overview.md) | End-to-end workflow, stage by stage, with inputs and outputs |
| [docs/02_directory_layout.md](docs/02_directory_layout.md) | What lives where, naming conventions, environment variables |
| [docs/03_file_mapping.md](docs/03_file_mapping.md) | Every script: original path → new path, and what changed |
| [docs/04_script_inputs_outputs.md](docs/04_script_inputs_outputs.md) | Inputs, outputs and the section-by-section figure map |
| [docs/05_software_and_environment.md](docs/05_software_and_environment.md) | Software versions, conda environments, R packages, databases |

## Reproducing the analyses

**1. Upstream pipelines.** The scripts under `pipelines/` are the job scripts used for the
production runs (`#SBATCH` headers included). They are reproduced as executed, except that the
site-specific absolute paths were replaced by variables: each script starts with a
`# Site-specific settings` block listing them (for example `ORF_DIR`, `SARG_DB`, `DIAMOND_DIR`,
`IN_DIR`). Set those variables and submit in the order given in
[docs/01](docs/01_pipeline_overview.md).

**2. R analyses.** The scripts under `analyses/` require R >= 4.2 and the packages listed in
[docs/05](docs/05_software_and_environment.md). Input tables are resolved through `CS_DATA`
and figures are written to `CS_WORK`; both default to the repository layout and can be
overridden without editing any script:

```r
Sys.setenv(CS_DATA = "/path/to/data", CS_WORK = "/path/to/results")
source("analyses/02_metagenome_resistome/01_resistome_analysis.R")
```

Run a script from top to bottom: the sections inside one notebook reuse the objects created
above them (for example the ARG type table loaded in section 3 is used again by the NMDS, the
heat maps and the risk analysis).

## Data availability

The small tables needed to follow the analyses (ASV table, SINTAX taxonomy, MAG abundance,
GTDB-Tk summary, plasmid counts, sample metadata, LEfSe input tables) are included under
[`data/`](data/README.md). The larger intermediate tables — ARGs-OAP output, diamond hits and
the MAG × gene matrices — are not distributed; they are regenerated by the scripts in
`pipelines/` from the raw reads. The raw sequence data are available from the corresponding
author on request, as stated in the manuscript. All databases used (SARG, BacMet, AromaDeg, MGE,
MiDAS 4.8.1, GTDB R214) are publicly available.

## Citation

```bibtex
@article{Li2025coselection,
  title   = {Understanding pollutant-driven shifts of antibiotic resistome in activated
             sludge: A lab-scale study},
  author  = {Li, Pengwei and Yin, Zirui and Ye, Lin},
  journal = {Journal of Hazardous Materials},
  volume  = {487},
  pages   = {137199},
  year    = {2025},
  doi     = {10.1016/j.jhazmat.2025.137199}
}
```

If you use the scripts, please cite the article above. The repository is also referenced in the
manuscript (Methods 2.3).
