#!/bin/bash -x
#SBATCH --mail-user=bougous
#SBATCH --mail-type=END
#SBATCH --output=slurm-AEG4844-kmer_count-k51-b1500-hdd-%j.out 
#SBATCH -N 1
#SBATCH -p workq
#SBATCH -t 12:00:00
#SBATCH -A k1616
#SBATCH -J AEG4844-kmer_count-k51-b1500-hdd

acc="AEG4844"
bins="1500"
outpath="output/debug/hdd"

mkdir -p ${outpath}

#####################
#echo "start decomp: "`date`
#pigz -vdf -c  /project/k1616/Guotai/gwas/rawdata/sequence_data/AEG4844/AEG4844_1.fq.gz > data/AEG4844_1.fq &
#pigz -vdf -c  /project/k1616/Guotai/gwas/rawdata/sequence_data/AEG4844/AEG4844_2.fq.gz > data/AEG4844_2.fq &
#wait
#echo "end decomp: "`date`
#####################

#####################
srun ./build/kmer_count_k51_c200k_debug $acc $bins $outpath &
#./build/kmer_count_k51_c200k $acc $bins $outpath &
wait
sleep 60

echo "kmer counter DONE! "`date`
#####################



