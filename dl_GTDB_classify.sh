#!/bin/bash
#SBATCH -o job.%j.%N.out 
#SBATCH --partition=fat

#SBATCH -J GTDB_lpw
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mail-type=end

source activate gtdbtk2.3

binning_path=/lustre/home/liutang_faculty/83LPW/02RDN/12drep/res_merged/dereplicated_genomes
#out_path=/lustre/home/liutang_faculty/11basic_lake/02DB/07GTDB



start_time=$(date +%s)
 
mkdir res
gtdbtk classify_wf --cpus $SLURM_NTASKS --pplacer_cpus  $SLURM_NTASKS  -x fa \
    --genome_dir ${binning_path}/  \
    --mash_db /lustre/home/liutang_faculty/02database/gtdbtk_release214/mash  \
    --out_dir res

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "running time is $(($cost_time/60)) min" >> run_time_gtdb_${SLURM_JOBID}.log



#mkdir ${out_path}/res

#gtdbtk classify_wf --cpus 14 \
#--genome_dir ${binning_path}/res_merged/dereplicated_genomes \
#--out_dir ${out_path}/res \
#-x fa \
#--pplacer_cpus 1

