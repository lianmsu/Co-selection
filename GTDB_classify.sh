#!/bin/bash

source activate gtdbtk2.3
binning_path=/lustre/home/liutang_faculty/83LPW/02RDN/12drep/res_merged/dereplicated_genomes
start_time=$(date +%s)
 
mkdir res
gtdbtk classify_wf --cpus 28 --pplacer_cpus  28  -x fa \
    --genome_dir ${binning_path}/  \
    --mash_db ~/gtdbtk_release214/mash  \
    --out_dir res

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "running time is $(($cost_time/60)) min" >> run_time_gtdb_${SLURM_JOBID}.log


