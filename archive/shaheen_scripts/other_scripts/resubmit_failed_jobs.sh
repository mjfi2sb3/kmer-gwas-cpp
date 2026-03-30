#!/bin/bash
#
#
timestamp=`date +%s`

jobid=$1

if [ ! -f "$jobid" ]; then
        echo "Please provide a file containing job ids with job_ids_<timestep>.txt"
        exit 1
fi

for i in `awk '{print $2}' ${jobid} | xargs -i sacct -j {} -o "jobid,state" --noheader|& grep -v "COMPLETED"|awk '{print $1}'| sed 's#\.0##g'|sort|uniq`; do grep $i ${jobid} | awk '{print $1}'; done > accessions_failed_${timestamp}.txt

echo "#########################################################"
echo "# Total num of accessions: "`wc -l ${jobid}`
echo "# Temporary accession list: accessions_failed_${timestamp}.txt"
echo "# Number of failed accessions: "`wc -l accessions_failed_${timestamp}.txt`
echo "#########################################################"

printf "\n\n# resubmitting ...  \n\n"
 bash ./kmer_gwas_shaheen2_bb_wrapper.sh accessions_failed_${timestamp}.txt
 
