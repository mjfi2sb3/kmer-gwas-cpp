#!/bin/bash
#SBATCH --time=10:00
#SBATCH -A k1616
#SBATCH -N 1
#SBATCH -p debug

#DW persistentdw name=salim_test
#DW stage_out type=file destination=/scratch/bougous/kmer_gwas/output/ source=$DW_PERSISTENT_STRIPED_salim_test/output/test_bb_i3.tar.gz

cd $DW_PERSISTENT_STRIPED_salim_test/output/

tar zcvf test_bb_i3.tar.gz test_bb_i3


exit 0;
