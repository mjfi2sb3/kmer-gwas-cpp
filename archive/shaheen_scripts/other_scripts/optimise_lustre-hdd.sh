#!/bin/bash -x
#
#
p="AEG0397";
suff="_luste-hdd"
prevID="28529603"

for bin in 1000 1500 2000; do

        cd ./data;
        bash ./mksymlink.sh $p $bin $suff
        cd ..;

        f=${p}"_k51_c200k_bin"$bin$suff

        #currentID=`sbatch --mem=0 --kill-on-invalid-dep=yes --dependency=afterany:${prevID} --parsable -J ${f}_kmer_gwas --output=slurm-${f}-%j-fscratch.err ./kmer_gwas_shaheen2_lustre-hdd.sh $f ${bin}`
	echo sbatch --mem=0  -J ${f}_kmer_gwas --output=slurm-${f}-%j-fscratch.err ./kmer_gwas_shaheen2_lustre-hdd.sh $f ${bin}
        #prevID=${currentID}

done

# limit RAM to 120 like Shaheen
#f=${p}"_k51_c200k_bin500"
#sbatch --mem=120g --kill-on-invalid-dep=yes --exclusive --dependency=afterany:${prevID} --parsable -J ${f}_kmer_gwas --output=slurm-${f}-%j-fscratch.err ./kmer_gwas_ibex.sh $f 500
