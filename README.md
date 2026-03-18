## how to run locally
nextflow run main.nf -profile standard \
      --accessions_file accessions.txt \
      --data_dir ./data \
      --num_bins 5
      
## how to run using slurm: ibex or shaheen
nextflow run main.nf -profile slurm \
--accessions_file accessions.txt \
--num_bins 10 --count n --delimiter 'none'
