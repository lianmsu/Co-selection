#!/bin/bash
#SBATCH -o job.%j.%N.out
#SBATCH -J gtdbtk
#SBATCH --partition=fat
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.gtdbtk.%j.%N.out
# =============================================================================
# Co-selection -- taxonomic classification of the dereplicated MAGs (GTDB-Tk)
# =============================================================================
# Runs the classify_wf workflow (ANI screening + pplacer placement) on the
# dereplicated genome set, which also produces the tree used for the iTOL
# annotation in analyses/02_metagenome_resistome.
#
# Manuscript      : Methods 2.3 (MAG recovery and taxonomy)
# Main outputs    : res/ (GTDB-Tk summary tables and the reference tree)
# Original script : GTDB_classify.sh
# -----------------------------------------------------------------------------
# Site-specific settings -- export them (or edit the defaults) before running:
#
#   BIN_DIR        directory with the dereplicated genomes (*.fa)   [required]
#   GTDB_MASH_DB   GTDB-Tk mash database directory                  [required]
#   OUT_DIR        output directory                                 [res]
#   CPUS           CPUs for classify_wf and pplacer                 [28]
# =============================================================================
BIN_DIR="${BIN_DIR:?set BIN_DIR to the directory with the dereplicated genomes}"
GTDB_MASH_DB="${GTDB_MASH_DB:?set GTDB_MASH_DB to the GTDB-Tk mash database}"
OUT_DIR="${OUT_DIR:-res}"
CPUS="${CPUS:-28}"

source activate gtdbtk2.3

start_time=$(date +%s)

mkdir -p "${OUT_DIR}"
gtdbtk classify_wf --cpus "${CPUS}" --pplacer_cpus "${CPUS}" -x fa \
    --genome_dir "${BIN_DIR}/" \
    --mash_db "${GTDB_MASH_DB}" \
    --out_dir "${OUT_DIR}"

end_time=$(date +%s)
cost_time=$((end_time - start_time))
echo "running time is $((cost_time / 60)) min" >> "run_time_gtdb_${SLURM_JOBID}.log"
