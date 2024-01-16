#!/bin/bash
#---------------Script SBATCH - NLHPC ----------------
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

source $HOME/miniconda3/bin/activate
source activate bbmap

MITOREF=/home/fleon/tortu-genomes/mitogenome_03/
TRIMDIR=/home/fleon/tortu-genomes/trimommatic_02
SPLITDIR=/home/fleon/tortu-genomes/mitogenome_03

list="GTFFS21_101_CKDL230023566-1A_H75N5DSX7_L4 GTFFS21_134_CKDL230023566-1A_H75N5DSX7_L4 GTFFS21_136_CKDL230023566-1A_H75N5DSX7_L4 GTFFS21_171_CKDL230023566-1A_H75N5DSX7_L4 GTFFS21_281_CKDL230023566-1A_H75N5DSX7_L4 GTFFS21_282_CKDL23$

for SAMPLE in $list
do


bbsplit.sh ref=$MITOREF/mitogenome.fasta in1=$TRIMDIR/${SAMPLE}_.trimmed_1P.fq in2=$TRIMDIR/${SAMPLE}_.trimmed_2P.fq basename=${SAMPLE}_mt_%.fq out1=$SPLITDIR/${SAMPLE}_nomt_R1.fq out2=$SPLITDIR/${SAMPLE}_nomt_R2.fq

done
