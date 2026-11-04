# Directory layout

```
Co-selection/
├── README.md / README_zh.md
├── pipelines/                       # upstream bioinformatics (shell + one python helper)
│   ├── 01_16S_amplicon/
│   │   └── 01_usearch_16S_amplicon.sh
│   ├── 02_metagenome_binning/
│   │   └── 02_gtdbtk_classify.sh
│   ├── 03_resistance_gene_annotation/
│   │   ├── 03_diamond_blastp_annotation.sh
│   │   ├── 04_diamond_results_to_tables.sh
│   │   ├── 05_gene_counts_per_genome.sh
│   │   └── turn_to_table.py
│   └── 04_gene_abundance/
│       └── 06_args_oap_gene_abundance.sh
├── analyses/                        # downstream statistics and figures (R)
│   ├── 01_16S_community/
│   │   └── 01_16S_community_analysis.R
│   ├── 02_metagenome_resistome/
│   │   └── 01_resistome_analysis.R
│   └── 03_MAG_gene_annotation/
│       └── 01_MAG_gene_cooccurrence.R
├── data/                            # small input tables, see data/README.md
└── docs/
```

## Conventions

* **Two layers.** `pipelines/` holds everything that produces a table from
  sequencing data (Slurm job scripts and one python helper); `analyses/` holds
  everything that turns a table into a statistic or a figure (R).
* **Numbering.** `pipelines/` is numbered continuously in execution order
  (`01` … `06`, plus the unused `turn_to_table.py` helper). `analyses/` is
  numbered per module and each module lists the figures it produces in
  `docs/01_pipeline_overview.md`.
* **File headers.** Every script starts with a header stating its purpose, the
  manuscript item it produces, its main outputs and its path in the original
  working copy. The two long R scripts additionally carry a section index with
  the line numbers of every analysis block.
* **No absolute path.** The production runs used absolute paths: a personal
  Windows data folder for the R notebooks and home/scratch directories on the
  cluster for the shell scripts. All of them were replaced by variables, so no
  user name, account, machine or personal directory appears in the repository.

## Environment variables

| Variable | Used by | Default | Meaning |
| --- | --- | --- | --- |
| `CS_HOME` | R scripts | `.` | repository root |
| `CS_DATA` | R scripts | `$CS_HOME/data` | root of the input tables |
| `CS_WORK` | R scripts | `$CS_HOME/results` | where figures/tables are written (created if missing) |
| `BIN_DIR` | `02_gtdbtk_classify.sh` | — | dereplicated genomes |
| `GTDB_MASH_DB` | `02_gtdbtk_classify.sh` | — | GTDB-Tk mash database |
| `ORF_DIR`, `OUT_DIR` | `03_diamond_blastp_annotation.sh` | — | Prodigal ORFs, diamond output |
| `SARG_DB`, `MRG_DB`, `MGE_DB`, `AROMADEG_DB` | `03_diamond_blastp_annotation.sh` | — | diamond databases |
| `DIAMOND_DIR`, `OUT_ROOT` | `04_diamond_results_to_tables.sh`, `05_gene_counts_per_genome.sh` | — | diamond output, merged tables |
| `SARG_STRUCT`, `MRG_STRUCT`, `ADG_STRUCT` | `04_diamond_results_to_tables.sh` | — | database structure files |
| `IN_DIR`, `OUT_ROOT`, `SARG_DB`, `ADG_DB`, `MRG_DB`, `MGE_DB` (+ `*_STRUCT`) | `06_args_oap_gene_abundance.sh` | — | clean reads, output root, databases |

Every script also reads `THREADS` / `CPUS`, `PYTHON` and the trim parameters of
the amplicon pipeline from the environment; the defaults are the values used in
the production runs.

## Where the outputs go

* R figures and summary tables: `results/` (git-ignored, created on the fly).
* LEfSe input tables and SPEC-OCCU output: `data/16S/` — the original notebook
  wrote them next to its input tables, and they are the files that the external
  LEfSe run consumes.
* `pipelines/` writes into the directories given by `OUT_DIR`, `OUT_ROOT` and
  `DIAMOND_DIR`, so nothing is written inside the repository by default.
