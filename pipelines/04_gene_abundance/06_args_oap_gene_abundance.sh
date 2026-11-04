#!/bin/bash
#SBATCH -o job.%j.%N.out
#SBATCH -J args_oap
#SBATCH --partition=fat
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
# =============================================================================
# Co-selection -- gene abundance quantification with ARGs-OAP v3.2.2
# =============================================================================
# Quantifies the abundance of ARGs, MRGs, ADGs and MGEs from the clean
# metagenomic reads (stage_one: alignment; stage_two: normalization to gene
# copies per cell).
#
# Manuscript      : Methods 2.3 (read-based annotation)
# Main outputs    : <OUT_ROOT>/<category>/normalized_cell.*.txt
# Original script : gene_abundance.sh
# -----------------------------------------------------------------------------
# Site-specific settings -- export them (or edit the defaults) before running:
#
#   IN_DIR        directory with the clean reads                   [required]
#   OUT_ROOT      root for the four result folders                 [required]
#   SARG_DB       SARG v3.2 short sub-database fasta               [required]
#   SARG_STRUCT   SARG v3.2 short sub-database structure file      [required]
#   ADG_DB        AromaDeg fasta                                   [required]
#   ADG_STRUCT    AromaDeg structure file                          [required]
#   MRG_DB        MRG fasta                                        [required]
#   MRG_STRUCT    MRG structure file                               [required]
#   MGE_DB        MGE fasta                                        [required]
#   MGE_STRUCT    MGE structure file                               [required]
#   THREADS       threads for both stages                          [28]
# =============================================================================
IN_DIR="${IN_DIR:?set IN_DIR to the directory with the clean reads}"
OUT_ROOT="${OUT_ROOT:?set OUT_ROOT to the output root}"
SARG_DB="${SARG_DB:?set SARG_DB to the SARG fasta}"
SARG_STRUCT="${SARG_STRUCT:?set SARG_STRUCT to the SARG structure file}"
ADG_DB="${ADG_DB:?set ADG_DB to the AromaDeg fasta}"
ADG_STRUCT="${ADG_STRUCT:?set ADG_STRUCT to the AromaDeg structure file}"
MRG_DB="${MRG_DB:?set MRG_DB to the MRG fasta}"
MRG_STRUCT="${MRG_STRUCT:?set MRG_STRUCT to the MRG structure file}"
MGE_DB="${MGE_DB:?set MGE_DB to the MGE fasta}"
MGE_STRUCT="${MGE_STRUCT:?set MGE_STRUCT to the MGE structure file}"
THREADS="${THREADS:-28}"

run_args_oap() {
    name="${1}"
    db="${2}"
    struct="${3}"
    out="${OUT_ROOT}/${name}"
    mkdir -p "${out}"
    echo "Step1 ${name}"
    args_oap stage_one -i "${IN_DIR}" -o "${out}" -t "${THREADS}" --database "${db}"
    echo "Step2 ${name}"
    args_oap stage_two -i "${out}" --database "${db}" -t "${THREADS}" --structure1 "${struct}"
}

source activate args_oap

run_args_oap ARG      "${SARG_DB}" "${SARG_STRUCT}"
run_args_oap ADG      "${ADG_DB}"  "${ADG_STRUCT}"
run_args_oap MRG      "${MRG_DB}"  "${MRG_STRUCT}"
run_args_oap MGE      "${MGE_DB}"  "${MGE_STRUCT}"
