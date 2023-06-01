#!/bin/bash
#SBATCH --mail-user=bougous
#SBATCH --mail-type=END
#SBATCH --output=slurm_out/AEG4844/slurm-AEG4844--kmer_count-k51-b1500-%j.out 
#SBATCH -N 1
#SBATCH -p workq
#SBATCH -t 8:00:00
#SBATCH -A k1616

#DW jobdw type=scratch access_mode=striped capacity=1TiB
#DW stage_in type=file source=/project/k1616/Guotai/gwas/rawdata/sequence_data/AEG4844/AEG4844_1.fq.gz destination=$DW_JOB_STRIPED/data/AEG4844_1.fq.gz
#DW stage_in type=file source=/project/k1616/Guotai/gwas/rawdata/sequence_data/AEG4844/AEG4844_2.fq.gz destination=$DW_JOB_STRIPED/data/AEG4844_2.fq.gz
#DW stage_in type=file source=/scratch/bougous/kmer_gwas/build/kmer_count_k51_c200k destination=$DW_JOB_STRIPED/kmer_count_k51_c200k
#DW stage_in type=file source=/scratch/bougous/kmer_gwas/kmer_gwas_shaheen2_bb.sh destination=$DW_JOB_STRIPED/kmer_gwas_shaheen2_bb.sh
#DW stage_out type=directory destination=/scratch/bougous/kmer_gwas/stageout_AEG4844 source=$DW_JOB_STRIPED/



cd $DW_JOB_STRIPED
chmod +x kmer_count_k51_c200k

acc="AEG4844"
bins="1500"
outpath="output"



mkdir -p ${outpath}

#####################
echo "start decomp: "`date`
pigz -fv  -d data/${acc}_1.fq.gz &
pigz -fv -d data/${acc}_2.fq.gz &
wait
echo "end decomp: "`date`
#####################

#####################
#srun ./kmer_count_k51_c200k $acc $bin $outpath &
./kmer_count_k51_c200k $acc $bins $outpath &
wait
sleep 60

echo "kmer counter DONE! "`date`
#####################



