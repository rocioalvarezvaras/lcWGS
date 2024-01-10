# lcWGS GREEN SEA TURTLE

---

# UPSTREAM ANALYSIS

# Step 1a: Quality Control of Raw reads - fastQC

FastQC performs a series of analyses, each represented by a module, to evaluate different aspects of data quality.
This tool takes as input one or more sequence files in various formats, such as FASTQ or BAM. In this case, we are going to use a fastq file.
- Input: fq.gz/fastq.gz files
- Output: .html and .zip files per sample (forward and reverse, separately)

Since we have numerous samples, we will use a loop that will help us analyze all the samples with the following script:

```bash
cat list.txt | while read a
do
  fastqc –o $GEN/${a} –t 1 –f fastq $GEN/${a}.gz
done
```
Details:
- The names of all the samples to be analyzed must be in the list.txt
- -o: the name that the output .html and .zip files will have
- -t: number of threads (cpu) to use to perform the analysis
- -f: to indicate that the input has fastq format, then you have provide the file path 

  
# Step 1b: Quality Control of Raw reads - multiQC
MultiQC is designed to aggregate and summarize results from multiple analysis tools, including FastQC, into a single, easy-to-read report.
To run multiqc we must be in the folder where our fastq files are located.
- Input: .fastqc files from all samples (forward and reverse, separately)
- Output: one .html file (all samples together)

Command:
``` bash
multiqc .

```
Details:
- By placing the '.' we are indicating that our input corresponds to all the fastqc files that are in our current directory

# Step 2: Trimming Raw reads - trimmomatic

- Input: fq.gz files
- Output: four files per sample, two unpaired (U) and two paired (P)
- Add list of adapters (Illumina)

_Execute script: step_2.sh_

## Run fastQC and multiQC post-trimming using paired files to re-evaluate the content of adapters

# Step 3: Separation of the mitogenome from the nuclear genome

- ref=$MITOREF: path to the reference mitogenome
- in1= path to reads (forward)
- in2= path to reads (reverse)
- basename= path and name of the file with the mitogenome (output name)
- out1= path and name of the nuclear reads (forward)
- out2= path and name of the nuclear reads (reverse)

_Execute script: step_3.sh_

Move the nuclear genome files to a new folder called "bbsplit_03" and the mitognomes to a new folder called "mitogen_03".

# Step 4: Genome indexing

Before aligning the reads to the reference genome, it is necessary to index the reference genome. 
The reference genome was used: GCF_015237465.2

-REF= path where the reference genome is
-p= Output database prefix [same as db file name]

_Execute script: step_4.sh_

The results should be in the “reference” folder and must be 5 files:
- cmydas_index.amb
- cmydas_index.ann
- cmydas_index.bwt
- cmydas_index.pac
- cmydas_index.sa
  
These files should be moved to a new folder called "bwa_04".

# Step 5: Genome alignment
