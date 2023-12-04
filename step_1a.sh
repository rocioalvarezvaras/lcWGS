#!/bin/bash
#---------------Script SBATCH - NLHPC ----------------
#SBATCH -J FastQC 
#SBATCH -p slims
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=2300
#SBATCH --mail-user=email
#SBATCH --mail-type=ALL
#SBATCH -t 24:2:5
#SBATCH -o FastQC_%j.out
#SBATCH -e FastQC_%j.err

GEN=/home/fleon/tortu-genomes 

source $HOME/miniconda3/bin/activate
source activate assembly

fastqc –o $GEN/ –t 1 –f fastq $GEN/*.gz
