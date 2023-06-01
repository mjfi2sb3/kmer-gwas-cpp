#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH -A k1616
#SBATCH --mem=0
#SBATCH -N 1
#SBATCH --mail-type=END
#SBATCH --mail-user=bougous


##acc="AEG0397_k51_c200k_bin500_lustre-hdd"; sbatch -J ${acc} --output=slurm-${acc}-%j.err kmer_gwas_shaheen2_lustre-hdd.sh ${acc} 1000

#module purge;
module load gcc/11.2.0

acc=$1
bin=$2
outpath=$3
#suff=$3

#cd data;
#bash ./mksymlink.sh $acc $bin  $suff
#cd ..;


#mkdir -p output

/usr/bin/time ./build/kmer_count_k51_c200k_v2 $acc $bin $outpath



