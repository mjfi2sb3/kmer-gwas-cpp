#!/bin/bash
#SBATCH --time=30:00
#SBATCH -A k1616
#SBATCH --mem=0
#SBATCH -N 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bougous
#SBATCH -p debug

#DW persistentdw name=salim_test
#DW stage_in type=directory source=/scratch/bougous/kmer_gwas destination=$DW_PERSISTENT_STRIPED_salim_test


##sbatch -J AEG0397_k51_c200k_b500_bb --output=slurm-AEG0397_k51_c200k_b500_bb-%j.err kmer_gwas_shaheen2.sh AEG0397 500

cd $DW_PERSISTENT_STRIPED_salim_test

chmod +x build/kmer_count_k51_c200k_v2

module load gcc/11.2.0

acc=$1
bin=$2
outpath=$3
#mkdir -p output

/usr/bin/time ./build/kmer_count_k51_c200k $acc $bin $outpath



