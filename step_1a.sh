#!/bin/bash
#SBATCH -J FastQC 
#SBATCH -p partition
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=2300
#SBATCH --mail-user=email
#SBATCH --mail-type=ALL
#SBATCH -o FastQC_%j.out
#-----Set relative path--------------------
GEN=/path/to/directory 
#-----Activate fastqc conda environment-----
source $HOME/miniconda3/bin/activate
activate assembly
#-----Command-------------------------------
fastqc –o $GEN/ –t 1 –f fastq $GEN/*.gz
