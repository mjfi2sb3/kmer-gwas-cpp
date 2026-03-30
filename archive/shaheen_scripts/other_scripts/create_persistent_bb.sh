#!/bin/bash
#SBATCH --time=30:00
#SBATCH -A k1616
#SBATCH --mem=0
#SBATCH -N 1
##SBATCH --mail-type=ALL
##SBATCH --mail-user=bougous
#SBATCH -p debug
#SBATCH -J create_persistent_space

#BB create_persistent name=salim_test capacity=300G access=striped type=scratch


exit 0
