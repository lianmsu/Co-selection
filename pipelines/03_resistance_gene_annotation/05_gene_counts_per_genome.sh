#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH -J gene_counts
#SBATCH --partition=fat
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# =============================================================================
# Co-selection -- per-genome counts of the four gene categories
# =============================================================================
# Counts the diamond hits of every genome and writes the matrix that is read by
# analyses/02_metagenome_resistome (section "per-MAG gene counts") and by the
# iTOL annotation step.
#
# Manuscript      : Fig. 3d, 3e (gene numbers per MAG)
# Main outputs    : <DIAMOND_DIR>/tongji2.txt
# Original script : diamond_res.sh (second half)
# -----------------------------------------------------------------------------
# Site-specific settings:
#
#   DIAMOND_DIR    output of 03_diamond_blastp_annotation.sh      [required]
# =============================================================================
DIAMOND_DIR="${DIAMOND_DIR:?set DIAMOND_DIR to the diamond output directory}"

summary="${DIAMOND_DIR}/tongji2.txt"
echo -e "genome\tARG\tMRG\tAromaDeg\tMGE" > "${summary}"

for folder in "${DIAMOND_DIR}"/*/
do
    [ -d "${folder}" ] || continue
    folder=$(basename "${folder}")
    echo "${folder}"

    ARG=$(awk 'END{print NR}' "${DIAMOND_DIR}/${folder}/"*ARG_out)
    MRG=$(awk 'END{print NR}' "${DIAMOND_DIR}/${folder}/"*MRG_out)
    AromaDeg=$(awk 'END{print NR}' "${DIAMOND_DIR}/${folder}/"*AromaDeg_out)
    MGE=$(awk 'END{print NR}' "${DIAMOND_DIR}/${folder}/"*MGE_out)

    echo -e "${folder}\t${ARG}\t${MRG}\t${AromaDeg}\t${MGE}" >> "${summary}"
done
