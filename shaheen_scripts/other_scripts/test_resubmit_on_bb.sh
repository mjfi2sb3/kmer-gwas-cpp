#!/bin/bash -x
#SBATCH -N 1
#SBATCH -p workq
#SBATCH -t 2:00
#SBATCH --mail-type=NONE
#SBATCH -A k1616

#DW jobdw type=scratch access_mode=striped capacity=1TiB
#DW stage_in type=file source=/project/k1616/Guotai/gwas/rawdata/sequence_data/AEG7977/AEG7977_1.fq.gz destination=$DW_JOB_STRIPED/data/AEG7977_1.fq.gz
#DW stage_in type=file source=/scratch/bougous/kmer_gwas/test_resubmit_on_bb.sh destination=$DW_JOB_STRIPED/test_resubmit_on_bb.sh
#DW stage_out type=directory destination=/scratch/bougous/kmer_gwas/test_resubmit_on_bb/ source=$DW_JOB_STRIPED/

cd /scratch/bougous/kmer_gwas/;
pwd;
###################################################
## resubmit if failed
sbatch --mem=262144 --kill-on-invalid-dep=yes --dependency=afternotok:$SLURM_JOBID -J re_$SLURM_JOB_NAME --mail-user=$USER --output=slurm-test_resubmit_on_bb-kmer_count-REQ-%j.out test_resubmit_on_bb.sh
sbatch --mem=262144 --kill-on-invalid-dep=yes --dependency=afternotok:$SLURM_JOBID -J re_$SLURM_JOB_NAME -A k1616 --mail-user=$USER --output=slurm-test_resubmit_on_bb-MEMHOG-REQ-%j.out --wrap="free -g; memhog 200G;"
###################################################


cd $DW_JOB_STRIPED


sleep 20;
exit 1




