#!/bin/bash

timestamp=`date +%s`

accessions=$1

if [ ! -f "$accessions" ]; then
	echo "Please provide a file containing accession names; one per line"
	exit 1
fi

# Prep some dirs
slurm_out="slurm_out/job_ids_${timestamp}"
mkdir -p ${slurm_out}

k=51 # kmer lenght is hardcoded. If must change, you will need to recompile code
bins=1500;
outpath="output/k${k}_b${bins}"

echo "#################################"
echo "Kmer length: ${k}"
echo "Number of Matrix Bins: ${bins}"
echo "Timestamp: ${timestamp} "
echo "Job ids are saved in:  job_ids_${timestamp}.log"
echo "Slurm output/error file location: "${slurm_out}
echo "#################################"

while read acc; do

	## create bb config
	sed "s#ACCESSION_REPLACE#${acc}#g" config/bb_tmpl.conf > config/bb_${acc}.conf

	## submit normal mem & full threading jobs
	id=`sbatch --parsable -J ${acc}_${k}_b${bins} --mail-user=$USER --output=${slurm_out}/slurm-${acc}-kmer_count-k${k}-b${bins}-%j.out --bbf=config/bb_${acc}.conf kmer_gwas_shaheen2_bb.sh ${acc} ${bins} ${outpath} ${slurm_out}`
	## submit himem & nomultihread jobs	
	#id=`sbatch --mem=256G --parsable -J ${acc}_${k}_b${bins} --mail-user=$USER --output=${slurm_out}/slurm-${acc}-kmer_count-k${k}-b${bins}-%j.out --bbf=config/bb_${acc}.conf kmer_gwas_shaheen2_bb_himem_nomultithread.sh ${acc} ${bins} ${outpath} ${slurm_out}`
	
	## record job ids	
	echo "$acc ${id}" >> job_ids_${timestamp}.log
done < $accessions
