#!/bin/bash
#SBATCH --time=30:00
#SBATCH -A k1616
#SBATCH --mem=0
#SBATCH -N 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bougous
#SBATCH -p debug

#BB create_persistent name=sb_test capacity=300G access=striped type=scratch
#DW jobdw type=persistent access_mode=striped capacity=1TiB
#DW stage_in type=directory source=/scratch/bougous/kmer_gwas destination=$DW_JOB_STRIPED
#DW stage_out type=directory destination=/scratch/bougous/kmer_gwas source=$DW_JOB_STRIPED/


##sbatch -J AEG0397_k51_c200k_b500_bb --output=slurm-AEG0397_k51_c200k_b500_bb-%j.err kmer_gwas_shaheen2.sh AEG0397 500

cd $DW_JOB_STRIPED
chmod +x build/kmer_count_k51_c200k_v2

acc=$1
bin=$2
outpath=$3

#mkdir -p output

/usr/bin/time ./build/kmer_count_k51_c200k $acc $bin output/k51_c200k_b${bin} $outpath

cd $outpath;
tar zcvf ${acc}.tar.gz ${acc} && rm -rf ${acc} 



