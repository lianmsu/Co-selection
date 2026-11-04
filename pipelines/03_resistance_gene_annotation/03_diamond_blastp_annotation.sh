#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH -J diamond_anno
#SBATCH --partition=fat
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
# =============================================================================
# Co-selection -- functional annotation of the MAG ORFs with diamond blastp
# =============================================================================
# Each genome folder produced by Prodigal is searched against the four
# databases used in the study: ARGs (SARG v3.2), MRGs (BacMet-based),
# ADGs (AromaDeg) and MGEs.
#
# Manuscript      : Methods 2.3 (ORF annotation with diamond v2.0.6)
# Main outputs    : <genome>/<genome>_{ARG,MRG,AromaDeg,MGE}_out
# Original script : diamond.sh
# -----------------------------------------------------------------------------
# Site-specific settings -- export them (or edit the defaults) before running:
#
#   ORF_DIR        Prodigal output, one folder per genome      [required]
#   OUT_DIR        diamond output root                          [required]
#   SARG_DB        SARG v3.2 diamond database                   [required]
#   MRG_DB         MRG diamond database                         [required]
#   MGE_DB         MGE diamond database                         [required]
#   AROMADEG_DB    AromaDeg diamond database                    [required]
#   THREADS        diamond threads per search                   [14]
# =============================================================================
ORF_DIR="${ORF_DIR:?set ORF_DIR to the Prodigal output directory}"
OUT_DIR="${OUT_DIR:?set OUT_DIR to the diamond output directory}"
SARG_DB="${SARG_DB:?set SARG_DB to the SARG diamond database}"
MRG_DB="${MRG_DB:?set MRG_DB to the MRG diamond database}"
MGE_DB="${MGE_DB:?set MGE_DB to the MGE diamond database}"
AROMADEG_DB="${AROMADEG_DB:?set AROMADEG_DB to the AromaDeg diamond database}"
THREADS="${THREADS:-14}"

cd "${ORF_DIR}"
for folder in *
do
    mkdir -p "${OUT_DIR}/${folder}"
    query="${folder}/protein_seq.fa"
    diamond blastp -k 1 -e 0.00001 -p "${THREADS}" -d "${SARG_DB}" \
        -q "${query}" -o "${OUT_DIR}/${folder}/${folder}_ARG_out" --id 80 --query-cover 80
    diamond blastp -k 1 -e 0.00001 -p "${THREADS}" -d "${MGE_DB}" \
        -q "${query}" -o "${OUT_DIR}/${folder}/${folder}_MGE_out" --id 80 --query-cover 80
    diamond blastp -k 1 -e 0.00001 -p "${THREADS}" -d "${MRG_DB}" \
        -q "${query}" -o "${OUT_DIR}/${folder}/${folder}_MRG_out" --id 80 --query-cover 80
    diamond blastp -k 1 -e 0.00001 -p "${THREADS}" -d "${AROMADEG_DB}" \
        -q "${query}" -o "${OUT_DIR}/${folder}/${folder}_AromaDeg_out" --id 80 --query-cover 80
done
