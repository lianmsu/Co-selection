# Pipeline overview

The study follows two sequencing arms that are analysed together: a 16S rRNA
gene amplicon arm (bacterial community) and a shotgun metagenomic arm (genes
and genomes). Both arms start from the same 15 activated-sludge samples.

## Experimental design

| Item | Value |
| --- | --- |
| Reactors | 21 lab-scale bioreactors (7 treatments × 3 replicates), 200 mL, MLSS 3 000 mg/L, pH 7, 25 °C |
| Treatments | blank control, Cr(VI), Cd(II), Cu(II), phenol, 4-aminophenol, 4-nitrophenol |
| Phases | Phase 1 (day 0–30) with pollutant, Phase 2 (day 30–60) without pollutant |
| Samples | 14 representative sludge samples (day 30 and day 60) + 1 seed sludge = 15 |
| Amplicon | 16S rRNA gene V4, Illumina MiSeq, 2 × 250 bp |
| Metagenome | Illumina HiSeq 4000, 2 × 150 bp |

Sample codes used everywhere in the scripts: `Seed`, `CG` (blank control),
`Cd`, `Cr`, `Cu` (heavy metals), `oh` (phenol), `PAP` (4-aminophenol) and
`PNP` (4-nitrophenol); the trailing `1`/`2` marks day 30 / day 60
(`Cd11` = Cd, day 30; `Cd21` = Cd, day 60).

## 16S rRNA gene arm

| Stage | Script | Input | Output |
| --- | --- | --- | --- |
| 1. Merge, trim, denoise, ASV table, rarefaction, SINTAX | `pipelines/01_16S_amplicon/01_usearch_16S_amplicon.sh` | paired-end fastq | `otutab.txt` (rarefied to 50 000 reads), `otus.fa`, `sintax.txt` |
| 2. Community composition, diversity, differential taxa, NCM, SPEC-OCCU | `analyses/01_16S_community/01_16S_community_analysis.R` | `data/16S/*`, `data/metadata.txt` | Fig. 6, Figs. S7–S8 |

## Metagenomic arm

| Stage | Script | Input | Output |
| --- | --- | --- | --- |
| 1. Read QC, assembly, binning (metaWRAP) and dereplication (dRep) | not included — standard metaWRAP/dRep commands, see `docs/05` | raw reads | 257 dereplicated MAGs |
| 2. Taxonomy and reference tree | `pipelines/02_metagenome_binning/02_gtdbtk_classify.sh` | dereplicated genomes | GTDB-Tk summary, tree for iTOL |
| 3. ORF prediction and diamond annotation | `pipelines/03_resistance_gene_annotation/03_diamond_blastp_annotation.sh` | Prodigal ORFs | `*_ARG_out`, `*_MRG_out`, `*_AromaDeg_out`, `*_MGE_out` |
| 4. Merge hits into genome × gene tables | `pipelines/03_resistance_gene_annotation/04_diamond_results_to_tables.sh`, `05_gene_counts_per_genome.sh` | diamond output + structure files | `*_gene_out1_table.csv`, `tongji2.txt` |
| 5. Read-based gene quantification (ARGs-OAP v3.2.2) | `pipelines/04_gene_abundance/06_args_oap_gene_abundance.sh` | clean reads | `normalized_cell.*.txt` for ARG/MRG/ADG/MGE |
| 6. Resistome analyses | `analyses/02_metagenome_resistome/01_resistome_analysis.R` | stages 2–5 | Fig. 1, Fig. 2, Fig. 3c, Fig. 5, Figs. S2–S6 |
| 7. Gene co-occurrence on MAGs | `analyses/03_MAG_gene_annotation/01_MAG_gene_cooccurrence.R` | stage 4 tables | Fig. 3d, Fig. 4, Fig. S4 |

## How the pieces connect

```
16S reads ─► usearch ─► ASV table + SINTAX taxonomy ─┐
                                                      ├─► analyses/01 ─► Fig. 6, S7, S8
metadata.txt ─────────────────────────────────────────┘

metagenomic reads ─► metaWRAP + dRep ─► MAGs
                                        ├─► GTDB-Tk ─────► taxonomy, tree ───────────┐
                                        ├─► Prodigal ─► diamond ─► MAG × gene tables ─┼─► analyses/02 ─► Fig. 1, 2, 3c, 5, S2-S6
                                        └─► ARGs-OAP ─► normalized abundance tables ──┘

MAG × gene tables ─► analyses/03 ─► Fig. 3d, 4, S4
```

The three R scripts are sequential notebooks: each one starts from the input
tables on disk, so they can be run independently, but the sections inside a
script reuse the objects created by the sections above them.

## Figure → script map

| Figure | Content | Script |
| --- | --- | --- |
| Fig. 1 | ARG abundance and diversity, ARG-type composition, NMDS | `analyses/02_metagenome_resistome/01_resistome_analysis.R` |
| Fig. 2 | Top-100 ARG heat map, shared/unique ARG subtypes | `analyses/02_metagenome_resistome/01_resistome_analysis.R` |
| Fig. 3a, 3b | ARG : ADG and ARG : MRG regressions | `analyses/02_metagenome_resistome/01_resistome_analysis.R` |
| Fig. 3c | Mantel tests (ASV, ARG, MRG, ADG) | `analyses/02_metagenome_resistome/01_resistome_analysis.R` |
| Fig. 3d, 3e | MAG tree with gene counts, Gammaproteobacteria subtree | `pipelines/02_*` + iTOL section of `analyses/02_*` |
| Fig. 4 | ARG–MRG and ARG–ADG networks, resistance mechanisms | `analyses/03_MAG_gene_annotation/01_MAG_gene_cooccurrence.R` |
| Fig. 5 | Risk levels, MGE Venn, ARG : MGE regression, plasmid content | `analyses/02_metagenome_resistome/01_resistome_analysis.R` |
| Fig. 6 | Richness, NMDS, top-50 genus heat map | `analyses/01_16S_community/01_16S_community_analysis.R` |
| Fig. S2–S3 | ARG types, resistance mechanisms | `analyses/02_metagenome_resistome/01_resistome_analysis.R` |
| Fig. S4 | Spearman correlation of MAG-level genes | `analyses/03_MAG_gene_annotation/01_MAG_gene_cooccurrence.R` |
| Fig. S5–S6 | MRG and ADG composition | `analyses/02_metagenome_resistome/01_resistome_analysis.R` |
| Fig. S7 | LEfSe and SPEC-OCCU | `analyses/01_16S_community/01_16S_community_analysis.R` |
| Fig. S8 | Neutral community model and C-score | `analyses/01_16S_community/01_16S_community_analysis.R` |
