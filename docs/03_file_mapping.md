# File mapping (original repository → this repository)

The previous repository was flat: three R notebooks and five shell scripts at
the top level, with the input tables stored inside the two analysis folders.
Nothing was renamed semantically, so the mapping below is mostly a move.

## Scripts

| Original path | New path | Change |
| --- | --- | --- |
| `16S_R_analysis/16s_bishe.r` | `analyses/01_16S_community/01_16S_community_analysis.R` | path variables, header + section index |
| `metagenomic_R_analysis/meta_bishe.R` | `analyses/02_metagenome_resistome/01_resistome_analysis.R` | path variables, header + section index, two `ggsave` calls now write to `CS_WORK` |
| `metagenomic_R_analysis/maggene.R` | `analyses/03_MAG_gene_annotation/01_MAG_gene_cooccurrence.R` | path variables, header + section index |
| `16SrRNA.sh` | `pipelines/01_16S_amplicon/01_usearch_16S_amplicon.sh` | settings block, shebang added, `usearch11` typo fixed, trim/rarefaction parameters exposed |
| `GTDB_classify.sh` | `pipelines/02_metagenome_binning/02_gtdbtk_classify.sh` | settings block, output directory becomes a variable |
| `diamond.sh` | `pipelines/03_resistance_gene_annotation/03_diamond_blastp_annotation.sh` | settings block, four database paths and thread count become variables |
| `diamond_res.sh` | `pipelines/03_resistance_gene_annotation/04_diamond_results_to_tables.sh` **and** `05_gene_counts_per_genome.sh` | split: the merging/annotation half and the per-genome count half; the three repeated ARG/MRG/ADG blocks became one loop; `turn_to_table.py` is now called through its own directory |
| `turn_to_table.py` | `pipelines/03_resistance_gene_annotation/turn_to_table.py` | portable shebang, `os.path.basename` instead of splitting on `/`, `main()` wrapper |
| `gene_abundance.sh` | `pipelines/04_gene_abundance/06_args_oap_gene_abundance.sh` | settings block, the four database runs became one function |

## Input tables

| Original path | New path |
| --- | --- |
| `16S_R_analysis/oturaw.txt` | `data/16S/oturaw.txt` |
| `16S_R_analysis/sintax.txt` | `data/16S/sintax.txt` |
| `16S_R_analysis/asv_lefse_category.txt` | `data/16S/asv_lefse_category.txt` |
| `16S_R_analysis/asv_lefse_30d.txt` | `data/16S/asv_lefse_30d.txt` |
| `16S_R_analysis/asv_lefse_60d.txt` | `data/16S/asv_lefse_60d.txt` |
| `16S_R_analysis/kw.list.txt` | `data/16S/kw.list.txt` |
| `16S_R_analysis/spec_occu_specialist_all.txt` | `data/16S/spec_occu_specialist_all.txt` |
| `16S_R_analysis/metadata.txt` | `data/metadata.txt` (shared by the two arms) |
| `metagenomic_R_analysis/metadata2.txt` | `data/metadata2.txt` |
| `metagenomic_R_analysis/bins_abundance_tpm.tsv` | `data/MAGs/bins_abundance_tpm.tsv` |
| `metagenomic_R_analysis/gtdbtk.bac120.summary.tsv` | `data/MAGs/gtdbtk.bac120.summary.tsv` |
| `metagenomic_R_analysis/plasmid_result.txt` | `data/MAGs/plasmid_result.txt` |

## Dropped

| File | Reason |
| --- | --- |
| `16S_R_analysis/spec_occu_specialist_30d.txt`, `..._60d.txt` | partial outputs of the same SPEC-OCCU section; `spec_occu_specialist_all.txt` (the file the script writes) is kept |
| `README.md` (original) | replaced by a full description of the study, the layout and the figure map |
| working title "Co-selection mechanism drives elevated antibiotic resistance in activated sludge: a laboratory-scale study" | the repository now uses the published title and DOI |

## What changed inside the scripts

Only paths and output locations were touched; the analysis code itself is
unchanged.

* **R.** Every `setwd('<data root>/...')` became
  `setwd(file.path(CS_DATA, "<subfolder>"))`, and the absolute paths that were
  used directly in `read.delim()`/`read.csv()` became
  `file.path(CS_DATA, ...)`. Figure and table writes that used a bare file name
  now go to `CS_WORK`. `data_folder <- ...` was removed (it was only used in
  commented-out code).
* **Shell.** The absolute home and scratch paths were replaced by the variables
  listed in `docs/02_directory_layout.md`, with the production values as defaults;
  the `#SBATCH` headers were kept.
* **Filenames.** One Chinese figure name (`随机森林丰度组合图.pdf`) became
  `random_forest_top_taxa.pdf`; the other output names are unchanged.

The header of `analyses/01_16S_community/01_16S_community_analysis.R` and
`analyses/02_metagenome_resistome/01_resistome_analysis.R` contains a section
index, because the original files were single sequential notebooks covering
several figures each; `docs/04_script_inputs_outputs.md` expands the same
information into a table.
