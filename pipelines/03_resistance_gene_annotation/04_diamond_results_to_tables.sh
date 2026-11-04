#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH -J diamond_res
#SBATCH --partition=fat
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# =============================================================================
# Co-selection -- merge the diamond hits into per-gene tables
# =============================================================================
# For ARGs, MRGs and ADGs the per-genome diamond output is concatenated into a
# single table, joined to the database structure file (gene -> type / subtype)
# and converted into a genome x gene count table by turn_to_table.py.
#
# Manuscript      : Methods 2.3 (gene presence on MAGs)
# Main outputs    : <out>/arg/arg_gene_out1, arg_gene_out_table.csv (and mrg, adg)
# Original script : diamond_res.sh
# -----------------------------------------------------------------------------
# Site-specific settings -- export them (or edit the defaults) before running:
#
#   DIAMOND_DIR    output of 03_diamond_blastp_annotation.sh      [required]
#   OUT_ROOT       root for the merged tables                     [required]
#   SARG_STRUCT    SARG v3.2 full structure file                  [required]
#   MRG_STRUCT     MRG structure file                             [required]
#   ADG_STRUCT     AromaDeg structure file                        [required]
#   PYTHON         python interpreter used by turn_to_table.py    [python3]
# =============================================================================
DIAMOND_DIR="${DIAMOND_DIR:?set DIAMOND_DIR to the diamond output directory}"
OUT_ROOT="${OUT_ROOT:?set OUT_ROOT to the output directory}"
SARG_STRUCT="${SARG_STRUCT:?set SARG_STRUCT to the SARG structure file}"
MRG_STRUCT="${MRG_STRUCT:?set MRG_STRUCT to the MRG structure file}"
ADG_STRUCT="${ADG_STRUCT:?set ADG_STRUCT to the AromaDeg structure file}"
PYTHON="${PYTHON:-python3}"

# 脚本所在目录，用于定位 turn_to_table.py
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

for gene in arg mrg adg
do
    case "${gene}" in
        arg) tag="ARG"      ; mapping="${SARG_STRUCT}" ;;
        mrg) tag="MRG"      ; mapping="${MRG_STRUCT}"  ;;
        adg) tag="AromaDeg" ; mapping="${ADG_STRUCT}"  ;;
    esac

    outpath="${OUT_ROOT}/${gene}"
    mkdir -p "${outpath}"

    # 合并所有基因组的命中，并在最后一列追加基因组名称
    merged="${outpath}/merged_${gene}_output_file"
    : > "${merged}"
    find "${DIAMOND_DIR}" -name "*${tag}_out" -type f | while read -r hit_file
    do
        genome=$(basename "$(dirname "${hit_file}")")
        awk -v foldername="${genome}" '{print $0 "\t" foldername}' "${hit_file}" >> "${merged}"
    done

    # 按 structure 文件补上基因的功能分类
    awk -F'\t' 'NR==FNR{a[$2]=$4; next} {print $0 "\t" (a[$2] ? a[$2] : "")}' \
        "${mapping}" "${merged}" > "${outpath}/${gene}_gene_out1"

    # 长表转宽表（genome x gene），输出 <gene>_gene_out1_table.csv
    ( cd "${outpath}" && "${PYTHON}" "${SCRIPT_DIR}/turn_to_table.py" "${gene}_gene_out1" )
done
