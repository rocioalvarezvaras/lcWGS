# lcWGS GREEN SEA TURTLE

-------------


# UPSTREAM ANALYSIS


# Step 1a: QC raw reads-fastQC


- Input: fq.gz files
- Output: .html and .zip files per sample (forward and reverse, separately) 



Script:

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

*sbatch script



# Step 1b: QC raw reads-multiQC

- Input: .html files all samples (forward and reverse, separately)
- Output: one .html file (all samples together)


Script:

multiqc .



# Step 2: Trimming raw reads-trimmomatic

- Input: fq.gz files
- Output: four files per sample, two unpaired (U) and two paired (P)
- Add list of adapters (Illumina)


Script:

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

GEN=/home/fleon/tortu-genomes
TRIM=/home/fleon/tortu-genomes/trimommatic_02
ADAPT=/home/fleon/Adapt

source $HOME/miniconda3/bin/activate
source activate assembly

cat list.txt | while read a
do
trimmomatic PE -threads 5 -phred33 $GEN/${a}1.fq.gz $GEN/${a}2.fq.gz -baseout $TRIM/${a}.trimmed.fq ILLUMINACLIP:$ADAPT/TruSeq2OV$
done

*nohup ./trim_loop.sh > trimming_rav.log 2>&1  </dev/null &


## Run fastQC and multiQC post-trimming using paired files to re-evaluate the content of adapters #












