#!/bin/bash
#SBATCH -o job.%j.%N.out
#SBATCH -J usearch16S
#SBATCH --partition=fat
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
# =============================================================================
# Co-selection -- 16S rRNA gene amplicon processing (usearch)
# =============================================================================
# Merge paired-end reads, trim primers, denoise with UNOISE3, build the ASV
# table, rarify to a fixed depth and classify the ASVs with SINTAX/MiDAS.
#
# Manuscript      : Methods 2.3 (16S rRNA gene sequencing)
# Main outputs    : merged.fq, otus.fa, otutab.txt, sintax.txt
# Original script : 16SrRNA.sh
# -----------------------------------------------------------------------------
# Site-specific settings -- export them (or edit the defaults) before running:
#
#   USEARCH        usearch v10 binary                       [./usearch]
#   SINTAX_DB      MiDAS 4.8.1 reference fasta              [required]
#   PREFIX_TRIM    bases removed from the 5' end            [27]
#   SUFFIX_TRIM    bases removed from the 3' end            [28]
#   MAXEE          expected-error threshold                 [1.0]
#   RAREFY_DEPTH   reads per sample for the rarefied table  [50000]
#   SINTAX_CUTOFF  SINTAX confidence threshold              [0.5]
#
# The script works on all *_R1.raw.fastq files in the current directory.
# =============================================================================
USEARCH="${USEARCH:-./usearch}"
SINTAX_DB="${SINTAX_DB:?set SINTAX_DB to the MiDAS 4.8.1 reference fasta}"
PREFIX_TRIM="${PREFIX_TRIM:-27}"
SUFFIX_TRIM="${SUFFIX_TRIM:-28}"
MAXEE="${MAXEE:-1.0}"
RAREFY_DEPTH="${RAREFY_DEPTH:-50000}"
SINTAX_CUTOFF="${SINTAX_CUTOFF:-0.5}"

# 16S Usearch
chmod +x "${USEARCH}"

# 双端拼接
"${USEARCH}" -fastq_mergepairs ./*_R1.raw.fastq -relabel @ -fastqout merged.fq

# 去引物 + 质量过滤
"${USEARCH}" -fastx_truncate merged.fq -stripleft "${PREFIX_TRIM}" \
    -stripright "${SUFFIX_TRIM}" -fastqout stripped.fq
"${USEARCH}" -fastq_filter stripped.fq -fastq_maxee "${MAXEE}" \
    -fastaout filtered.fa -relabel Filt
"${USEARCH}" -fastx_uniques filtered.fa -sizeout -relabel Uniq -fastaout uniques.fa

# UNOISE3 去噪得到 ASV
"${USEARCH}" -unoise3 uniques.fa -zotus asv.fa

head -n 4 asv.fa
sed 's/Zotu/ASV_/g' asv.fa > otusasv.fa
head -n 10 otusasv.fa
cp -f otusasv.fa otus.fa
head -n 10 otus.fa

# 生成 OTU/ASV 丰度表并按固定深度抽平
"${USEARCH}" -otutab merged.fq -db otus.fa -otutabout otutab_raw.txt

"${USEARCH}" -otutab_rare otutab_raw.txt -sample_size "${RAREFY_DEPTH}" \
    -output otutab.txt

# SINTAX 物种注释（MiDAS 4.8.1）
"${USEARCH}" -makeudb_usearch "${SINTAX_DB}" -output SINTAX_MiDAS4.8.1.udb
"${USEARCH}" -sintax otus.fa -db SINTAX_MiDAS4.8.1.udb -strand both \
    -tabbedout sintax.txt -sintax_cutoff "${SINTAX_CUTOFF}"
