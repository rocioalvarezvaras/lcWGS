#!/bin/bash
# --------------- Script SBATCH - NLHPC ----------------
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

# Activar entorno virtual
source $HOME/miniconda3/bin/activate
conda activate assembly

# Ejecutar FastQC
fastqc --outdir=$GEN --threads=1 --format=fastq $GEN/*.gz