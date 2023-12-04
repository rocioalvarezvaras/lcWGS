#!/bin/bash
#---------------Script SBATCH - NLHPC ----------------
#SBATCH -J BWA-ID
#SBATCH -p general
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=3000
#SBATCH -t 2:00:00
#SBATCH -o  Cmydas-index%j_%x.out
#SBATCH -e  Cmydas-index%j_%x.err

source $HOME/miniconda3/bin/activate
source activate assembly

REF=/home/fleon/tortu-genomes/reference

bwa index -p $REF/cmydas_index $REF/GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna
