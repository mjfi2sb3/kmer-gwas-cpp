#!/bin/bash

timestamp=`date +%s`
tmpdir=tmp_matrix_merge_${timestamp};

if (( $# < 3)); then
	echo ""
	echo "USAGE: `basename $0` <accession_list.txt> <kmer_count_output_path> <NUM_BINS> <MIN_KMER_OCCUR_ACROSS_ACCESSIONS=6> <BATCH_SIZE_PER_NODE=10>"
	echo ""
	exit 0;
fi



## batch_size = number of indices to run on a single node
acc_list=$1
if [ -z "$acc_list" ] || [ ! -f "$acc_list" ] ; then
	echo "Provide a valid file accession list; one accession per line"
	exit 1;
fi


## SET kmer_count_output_path
kmer_count_path=$2
if [ ! -d "$kmer_count_path" ] ; then
	echo "<<"$kmer_count_path">> kmer_count output dir doesn't exist!!"
	exit 1;
fi

## SET NUM BINS
num_bins=$3
if [ -z "$num_bins" ] || (( $num_bins < 1 )) ; then
	echo "Enter a valid number of bins [ NUM_BINS > 0 ]"
	exit 1;
fi

## SET MIN_OCCUR
min_occur=$4
if [ -z "$min_occur" ] ; then
	min_occur=6
fi


## SET BATCH_SIZE
batch_size=$5
if [ -z "$batch_size" ] ; then
	batch_size=10
fi

echo "########################################################"
printf "# LIST OF ACCESSIONS [FILE]:\t"${acc_list}"\n"
printf "# NUMBER OF ACCESSIONS:\t\t"`wc -l $acc_list`
printf "\n"

printf "# KMER_COUNT OUT PATH:\t\t"$kmer_count_path"\n"

printf "# NUM BINS PER ACCESSION:\t"$num_bins"\n"

printf "# MIN KMER OCCUR:\t\t"$min_occur"\n"

printf "# BATCH_SIZE PER NODE:\t\t"$batch_size"\n"

printf "# OUTPUT:\t\t\t./matrix\n"

printf "# SLURM OUT:\t\t\t$tmpdir\n"

echo "########################################################"

########################
########################
check_bins="n"
if [ "$check_bins" == "y" ]; then
	printf "\n\nChecking number of bins per accession in $kmer_count_path ..."
	while read acc; do
		num=`ls $kmer_count_path/${acc}/*_nr.tsv|wc -l`
		if [ ${num} -ne ${num_bins} ]; then
			printf "\n\tACC ${acc} has ${num}...NOT ${num_bins}";
			exit 1;
		fi
		printf "\n${acc} ... OK"
	done < $acc_list
	printf "\tDONE!\n"
fi
########################
########################
i=0
TOPINDEX=$(($num_bins - 1))

#timestamp=`date +%s`
#tmpdir=tmp_matrix_merge_${timestamp};
mkdir -p ${tmpdir};
while [[ $i -le ${TOPINDEX} ]]; do
	submitscript="${tmpdir}/submit_${i}.sh"
	echo "#!/bin/bash" > $submitscript
	chmod +x $submitscript	
	a=0
	while [[ $a -lt 10 ]]; do
		echo "/usr/bin/time -v ./build/matrix_merge_201acc ${kmer_count_path} ${acc_list} $(($a+$i)) ${min_occur} &" >> $submitscript
		((a = a+1))
	done
	echo "wait;" >> $submitscript

	## compress output
	a=0
	while [[ $a -lt 10 ]]; do
		echo "pigz -9 -p 3 matrix/"$(($a+$i))"_m.tsv &" >> $submitscript
		((a = a+1))
	done
        echo "wait;" >> $submitscript
	# submit job
	sbatch -Q -J merge_matrix_${i} -A k1616 --mail-type=END --mail-user=bougous@kaust.edu.sa -t 5:00:00 --output=${tmpdir}/slurm-${i}-%j.out ${submitscript}
	((i = i + $batch_size))
done
echo "JOBS SUBMITTED!!"
#echo "rm ${tmpdir}"

