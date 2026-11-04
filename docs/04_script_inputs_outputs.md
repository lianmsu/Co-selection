# Script inputs, outputs and section map

All R paths are relative to `CS_DATA` (`data/` by default); the pipeline paths
are given by the variables listed in `docs/02_directory_layout.md`.

## pipelines

### `pipelines/01_16S_amplicon/01_usearch_16S_amplicon.sh`

* **Input** — `*_R1.raw.fastq` (paired-end 16S reads), `SINTAX_DB`
  (MiDAS 4.8.1 fasta).
* **Output** — `merged.fq`, `stripped.fq`, `filtered.fa`, `uniques.fa`,
  `asv.fa` / `otus.fa`, `otutab_raw.txt`, `otutab.txt` (rarefied to
  `RAREFY_DEPTH`, 50 000), `sintax.txt`.
* **Parameters** — `-fastq_maxee 1.0`, `-stripleft 27 -stripright 28`,
  `-unoise3`, `-otutab_rare -sample_size 50000`, `-sintax_cutoff 0.5`.

### `pipelines/02_metagenome_binning/02_gtdbtk_classify.sh`

* **Input** — dereplicated genomes (`BIN_DIR`, `*.fa`), GTDB-Tk mash database.
* **Output** — `OUT_DIR` with the GTDB-Tk summary tables, the alignment and the
  tree used for the iTOL annotation; `run_time_gtdb_<jobid>.log`.

### `pipelines/03_resistance_gene_annotation/03_diamond_blastp_annotation.sh`

* **Input** — `ORF_DIR/<genome>/protein_seq.fa` (Prodigal), four diamond
  databases (`SARG_DB`, `MRG_DB`, `MGE_DB`, `AROMADEG_DB`).
* **Output** — `OUT_DIR/<genome>/<genome>_{ARG,MRG,AromaDeg,MGE}_out`.
* **Parameters** — `blastp -k 1 -e 1e-5 --id 80 --query-cover 80`.

### `pipelines/03_resistance_gene_annotation/04_diamond_results_to_tables.sh`

* **Input** — the diamond output above, `SARG_STRUCT`, `MRG_STRUCT`,
  `ADG_STRUCT`.
* **Output** — `OUT_ROOT/<gene>/merged_<gene>_output_file`,
  `OUT_ROOT/<gene>/<gene>_gene_out1`, `OUT_ROOT/<gene>/<gene>_gene_out1_table.csv`
  (`genome × gene` counts, written by `turn_to_table.py`).

### `pipelines/03_resistance_gene_annotation/05_gene_counts_per_genome.sh`

* **Input** — the diamond output above.
* **Output** — `DIAMOND_DIR/tongji2.txt`: genome, ARG, MRG, AromaDeg and MGE
  hit counts. This is the table behind the heat-map lanes of Fig. 3d/3e and it
  is read by `analyses/02` (`MAGs/tongji2.txt`).

### `pipelines/03_resistance_gene_annotation/turn_to_table.py`

* **Input** — long table whose last two columns are (gene, genome).
* **Output** — `<input_basename>_table.csv`, tab separated, genes in rows.

### `pipelines/04_gene_abundance/06_args_oap_gene_abundance.sh`

* **Input** — clean reads (`IN_DIR`), SARG/ADG/MRG/MGE fasta and structure
  files.
* **Output** — `OUT_ROOT/{ARG,ADG,MRG,MGE}/normalized_cell.*.txt` and the
  accompanying `rpkm.*`/`ppm.*` tables. `normalized_cell.*` is the
  gene-copies-per-cell normalization used in the manuscript.

## analyses

### `analyses/01_16S_community/01_16S_community_analysis.R` (1 163 lines)

**Input** — `data/16S/oturaw.txt`, `data/16S/sintax.txt`, `data/metadata.txt`.
**Output** — figures in `CS_WORK` (`random_forest_top_taxa.pdf`,
`random_forest_importance.csv`, `random_forest_kruskal_results.txt`), LEfSe
tables and `spec_occu_specialist_all.txt` in `data/16S/`.

| Line | Section | Manuscript |
| --- | --- | --- |
| 33 | Load ASV table, SINTAX taxonomy, rarefaction, alpha diversity, metadata | — |
| 95 | Phylum-level community composition (stacked bars) | Fig. 6e context |
| 143 | Class-level composition inside Pseudomonadota | Fig. 6e context |
| 196 | Gammaproteobacteria abundance (violin + Wilcoxon) | Fig. 6 |
| 263 | Alpha diversity, day 30 and day 60 (bar / box plots) | Fig. 6a, 6b |
| 357 | Heat map of the 50 most abundant genera | Fig. 6e |
| 399 | NMDS of the ASV community, ANOSIM and PERMANOVA | Fig. 6c, 6d |
| 445 | LEfSe input tables (`asv_lefse_*.txt`) | Fig. S7 |
| 493 | ASV co-occurrence network (igraph → Gephi) | — |
| 564 | Random forest, top-10 Kruskal-Wallis | — |
| 745 | Neutral community model (minpack.lm) | Fig. S8 |
| 837 | C-score null model of community assembly | Fig. S8 |
| 1081 | Specificity-occupancy (SPEC-OCCU), specialists/generalists | Fig. S7 |
| 1159 | Placeholders (RDA/VPA, SVM), unused | — |

Notes: `rrarefy()` is called without a seed, exactly as in the original
notebook, so the rarefied ASV table of a rerun differs by the rarefaction
draw. Set `set.seed()` before line 33 if bit-identical numbers are needed.

### `analyses/02_metagenome_resistome/01_resistome_analysis.R` (2 054 lines)

**Input** — `data/metadata.txt`, `data/metadata2.txt`,
`data/MAGs/{gtdbtk.bac120.summary.tsv,bins_abundance_tpm.tsv,tongji2.txt,genomeInformation.csv,plasmid_result.txt}`,
`data/{ARG,MRG,aromadeg,MGE}/normalized_cell.*.txt`, `data/ARG/{rpkm,ppm}.type.txt`,
`data/ARG/risk.txt`, `data/ARG/jili2/normalized_cell.Mechanism.group.txt`,
`data/reviewers/*`, `data/db/SARG/*`.
**Output** — figures in `CS_WORK` (`ARG_total_abundance.pdf`,
`p2_argtype_percentage.pdf`, …); the network section exports graph objects for
Gephi.

| Line | Section | Manuscript |
| --- | --- | --- |
| 46 | iTOL annotation files for the MAG tree | Fig. 3d |
| 194 | MAG taxonomy and abundance summary | — |
| 215 | ARG types (reads level) | Fig. 1b, Fig. S2 |
| 306 | ARG alpha diversity | Fig. 1b context |
| 382 | ARG NMDS, ANOSIM/PERMANOVA | Fig. 1c |
| 426 | MRG types (reads level) | Fig. S5 |
| 518 | ADG types | Fig. S6 |
| 617 | MGE types | Fig. 5b context |
| 719 | Total ARG abundance per treatment | Fig. 1a |
| 794 | Total MRG/ADG/MGE abundance | — |
| 878 | Shared and unique ARG subtypes (Venn) | Fig. 2b, 2c |
| 928 | Shared and unique MGE/MRG/ADG genes (Venn) | Fig. 5b |
| 975 | Per-MAG gene counts from diamond | Fig. 3d, 3e |
| 1006 | MAG abundance (TPM) | — |
| 1048 | MAG completeness, contamination, genome size | — |
| 1066 | Genome size of gene-carrying vs non-carrying MAGs | — |
| 1079 | Mantel tests (ASV, ARG, MRG, ADG) | Fig. 3c |
| 1228 | Gene abundance vs seed sludge | — |
| 1303 | ARG subtype heat map | — |
| 1347 | Multi-annotation pheatmap of ARGs | Fig. 2a |
| 1462 | Regression fits between gene categories | Fig. 3a, 3b, 5c |
| 1554 | ARG health-risk levels | Fig. 5a |
| 1614 | ARG resistance mechanisms | Fig. S3 |
| 1682 | ARG–MRG–ADG co-occurrence network | Fig. 4a, 4b |
| 1755 | Plasmid content of MAGs (PlasFlow) | Fig. 5d, 5e, 5f |
| 1890 | Reviewer additions: per-treatment ARG totals | — |

### `analyses/03_MAG_gene_annotation/01_MAG_gene_cooccurrence.R` (369 lines)

**Input** — `data/ARG/{sarg,mrg,aromadeg}_gene_out*table.csv` (from
`04_diamond_results_to_tables.sh`), `data/ARG/arg_mrg.csv`,
`data/ARG/arg_aroma.csv`, `data/ARG/normalized_cell.subtype.txt`,
`data/db/SARG/ARG_rank.txt`.
**Output** — co-occurrence graphs for Gephi/GraphML in `CS_WORK`, cumulative
abundance and mechanism plots.

| Line | Section | Manuscript |
| --- | --- | --- |
| 24 | Merge gene tables, Spearman matrix with BH correction | Fig. S4 |
| 63 | Distribution of the co-selected ARGs | — |
| 92 | ARG–MRG and ARG–ADG networks | Fig. 4a, 4b |
| 162 | Cumulative abundance of MRG-linked ARGs | Fig. 4c |
| 237 | Cumulative abundance of ADG-linked ARGs | Fig. 4d |
| 318 | Resistance mechanisms of MRG-linked ARGs (pie) | Fig. 4c |
| 342 | Resistance mechanisms of ADG-linked ARGs (pie) | Fig. 4d |
