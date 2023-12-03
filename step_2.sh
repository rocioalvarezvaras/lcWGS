#!/bin/bash
# --------------- Script SBATCH - NLHPC ----------------
#SBATCH -J Trimmomatic-rav
#SBATCH -p slims
#SBATCH -n 5
#SBATCH -c 1
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-user=ralvarez03@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 30:20:5
#SBATCH -o Trimmomatic_%j.out
#SBATCH -e Trimmomatic_%j.err

GEN=/home/fleon/tortu-genomes
TRIM=/home/fleon/tortu-genomes/trimmomatic_02
ADAPT=/home/fleon/Adapt

source $HOME/miniconda3/bin/activate
conda activate assembly

while read a; do
  trimmomatic PE -threads 5 -phred33 \
    ${GEN}/${a}1.fq.gz \
    ${GEN}/${a}2.fq.gz -baseout \
    ${TRIM}/${a}.trimmed.fq ILLUMINACLIP:${ADAPT}/TruSeq2OV
done < list.txt

# nohup ./trim_loop.sh > trimming_rav.log 2>&1 </dev/null &
