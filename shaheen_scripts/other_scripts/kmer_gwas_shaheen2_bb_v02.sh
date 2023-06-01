#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH -A k1616
#SBATCH -N 1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bougous

#DW jobdw type=scratch access_mode=striped capacity=1TiB
#DW stage_in type=directory source=/scratch/bougous/kmer_gwas/data destination=$DW_JOB_STRIPED/data
#DW stage_in type=directory source=/scratch/bougous/kmer_gwas/build destination=$DW_JOB_STRIPED/build
#DW stage_out type=directory destination=/scratch/bougous/kmer_gwas/output source=$DW_JOB_STRIPED/output/


##sbatch -J AEG0397_k51_c200k_b500_bb --output=slurm-AEG0397_k51_c200k_b500_bb-%j.err kmer_gwas_shaheen2.sh AEG0397 500

cd $DW_JOB_STRIPED

chmod +x build/kmer_count_k51_c200k

acc=$1
bin=$2
outpath=$3

mkdir -p output

/usr/bin/time ./build/kmer_count_k51_c200k $acc $bin $outpath &
wait
sleep 60

#echo "done deduping .... compressing now"

## v2
#pigz -p $SLURM_JOB_CPUS_PER_NODE $outpath/${acc}/*;
#tar cvf $outpath/${acc}.tar $outpath/${acc} && rm -rf $outpath/${acc};

## v3
#tar cvf $outpath/${acc}.tar $outpath/${acc} && rm -rf $outpath/${acc};

