#!/bin/bash -x
#SBATCH -N 1
#SBATCH -p workq
#SBATCH -t 8:00:00
#SBATCH --mail-type=END
#SBATCH -A k1616
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --hint=nomultithread


export OMP_NUM_THREADS=32


acc=$1
bins=$2
outpath=$3
slurm_out=$4


###################################################
## resubmit if failed
cd /scratch/bougous/kmer_gwas/;
sbatch --mem=256G --dependency=afternotok:$SLURM_JOBID -J re_$SLURM_JOB_NAME --mail-user=$USER --output=${slurm_out}/slurm-${acc}-kmer_count-b${bins}-REQ-%j.out --bbf=config/bb_${acc}.conf kmer_gwas_shaheen2_bb_himem_nomultithread.sh ${acc} ${bins} ${outpath} ${slurm_out}
###################################################

cd $DW_JOB_STRIPED
chmod +x kmer_count_k51_c200k

mkdir -p ${outpath}

#####################
echo "start decomp: "`date`
pigz -fv  -d data/${acc}_1.fq.gz &
pigz -fv -d data/${acc}_2.fq.gz &
wait
echo "end decomp: "`date`
#####################

#####################
srun --hint=nomultithread --cpus-per-task=${OMP_NUM_THREADS}  ./kmer_count_k51_c200k $acc $bins $outpath &
wait

#####################



