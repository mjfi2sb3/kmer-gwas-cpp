#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH -A k1616
#SBATCH --mem=0
#SBATCH -N 1
#SBATCH --mail-type=END
#SBATCH --mail-user=salim.bougouffa@kaust.edu.sa

#DW jobdw type=scratch access_mode=striped capacity=1TiB
#DW stage_in type=directory source=/scratch/bougous/kmer_gwas destination=$DW_JOB_STRIPED
#DW stage_out type=directory destination=/scratch/bougous/kmer_gwas source=$DW_JOB_STRIPED/

cd $DW_JOB_STRIPED
chmod +x build/accession_count_k51_200k_v2

module purge;
module load gcc/11.2.0

acc=$1
bin=$2

#pigz -d data/${acc}_1.fq.gz &
#pigz -d data/${acc}_2.fq.gz &
#wait

mkdir -p output

/usr/bin/time ./build/accession_count_k51_200k_v2 $acc $bin
