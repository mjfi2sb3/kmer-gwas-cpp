#!/bin/bash

timestamp=`date +%s`

max_jobs=200;
account="k1665"
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
num_lines=$(wc -l < "$accessions");
num_lines=$((num_lines*1))

echo "#################################"
echo "Kmer length: ${k}"
echo "Number of Matrix Bins: ${bins}"
echo "Number of Accessions: ${num_lines}"vvv
echo "Max Jobs Running: ${max_jobs}"
echo "Timestamp: ${timestamp} "
echo "Job ids are saved in:  job_ids_${timestamp}.log"
echo "Slurm output/error file location: "${slurm_out}
echo "running on login node: "`hostname`
echo "#################################"

for i in $(seq 1 $num_lines); do
	
	## create bb config
	acc=`sed -n "${i}p" ${accessions}`;

	while [[ $(squeue -u $USER --noheader | wc -l) -ge "$max_jobs" ]]; do
		sleep 10;
	done

	sed "s#ACCESSION_REPLACE#${acc}#g" config/bb_tmpl.conf > config/bb_${acc}.conf
	#printf "\n${i} [${acc}] <--> ";
	id=`sbatch --parsable -J ${acc}_${k}_b${bins} -N 1 -p workq  -t 8:00:00 --mail-type=FAIL --mail-user=$USER --output=${slurm_out}/slurm-${acc}-kmer_count-k${k}-b${bins}-%j.out -A ${account} --bbf=config/bb_${acc}.conf kmer_gwas_shaheen2_bb.sh ${acc} ${bins} ${outpath} slurm_out`;

	## record job ids	
	echo "$acc ${id}" >> job_ids_${timestamp}.log
	
done
exit;
#while read acc; do
#
#	## create bb config
#	sed "s#ACCESSION_REPLACE#${acc}#g" config/bb_tmpl.conf > config/bb_${acc}.conf
#
#	## submit job
#	echo "sbatch --parsable -J ${acc}_${k}_b${bins} -N 1 -p workq  -t 8:00:00 --mail-type=FAIL --mail-user=$USER --output=${slurm_out}/slurm-${acc}-kmer_count-k${k}-b${bins}-%j.out -A k1616 --bbf=config/bb_${acc}.conf kmer_gwas_shaheen2_bb.sh ${acc} ${bins} ${outpath} slurm_out"
#	#id=`sbatch --parsable -J ${acc}_${k}_b${bins} -N 1 -p workq  -t 8:00:00 --mail-type=FAIL --mail-user=$USER --output=${slurm_out}/slurm-${acc}-kmer_count-k${k}-b${bins}-%j.out -A k1616 --bbf=config/bb_${acc}.conf kmer_gwas_shaheen2_bb.sh ${acc} ${bins} ${outpath} slurm_out`
	
#	## record job ids	
#	#echo "$acc ${id}" >> job_ids_${timestamp}.log
#	prev_id=${id};
#	c=$((c+1))
#done < $accessions
