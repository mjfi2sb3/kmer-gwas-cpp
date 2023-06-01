#!/bin/bash
#SBATCH --mail-user=bougous
#SBATCH --mail-type=NONE
#SBATCH --output=slurm-test-32threads-%j.out 
#SBATCH -N 1
#SBATCH -p workq
#SBATCH -t 10:00
#SBATCH -A k1616
#SBATCH --hint=nomultithread


export OMP_NUM_THREADS=32

acc="AEG0397_5m"
bins="1500"
outpath="output/debug/"

#mkdir -p ${outpath}

srun exit 1


#####################
#sleep 15;
#####################

#####################
#srun ./kmer_count_k51_c200k $acc $bin $outpath &
#srun --hint=nomultithread --cpus-per-task=${OMP_NUM_THREADS}    ./build/kmer_count_k51_c200k $acc $bins $outpath &
#wait


#echo "kmer counter DONE! "`date`
#####################



