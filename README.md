# lcWGS Pacific Green Sea Turtle

---

# UPSTREAM ANALYSIS

# Step 1a: Quality control of raw reads - fastQC

FastQC performs a series of analyses, each represented by a module, to evaluate different aspects of data quality.
This tool takes as input one or more sequence files in various formats, such as FASTQ or BAM. In this case, we are going to use a fastq file.
- Input: fq.gz/fastq.gz files
- Output: .html and .zip files per sample (forward and reverse, separately)

Since we have numerous samples, we used a loop that will help us analyze all the samples with the following script:

Script:
```bash
#------Set relative path--------------------
RAW=/home/ralvarezv/raw_reads  #Directory where raw reads are located
#------Command------------------------------
cat list.txt | while read a
do
  fastqc –o $RAW/${a} –t 1 –f fastq $RAW/${a}.gz
done
```
Details:
- Names of all the samples (prefix) to be analyzed must be in the list.txt
- -o: the name that the output (.html and .zip files) will have
- -t: number of threads (cpu) to use to perform the analysis
- -f: to indicate that the input has fastq format, then you have provide the file path 

# Step 1b: Quality control of raw reads - multiQC

MultiQC was designed to aggregate and summarize results from multiple analysis tools (including FastQC) into a single, easy-to-read report.
To run multiqc we must be in the folder where our fastq files are located.
- Input: .fastqc files from all samples (forward and reverse, separately)
- Output: one .html file (all samples together)

Command:
``` bash
multiqc .
```
Details:
- By placing the '.' we are indicating that our input corresponds to all the fastqc files that are in our current directory.

# Step 2: Trimming raw reads - Trimmomatic

We used Trimmomatic to trim quality, remove adapter sequences, and filter out low-quality reads from raw sequencing data.
- Input: fq.gz files
- Output: four files per sample, two unpaired (U) and two paired (P)
- Add list of adapters (Illumina)

Script:
``` bash
#------Set relative path--------------------
RAW=                                  #Directory where raw data is located (input)
TRIM=/home/ralvarezv/trimmomatic_02/  #Directory where trimmed reads will be located (output)
ADAPT=/home/ralvarezv/adapters  #Directory where the adapters are located
#------Command------------------------------
cat list.txt | while read a
do
  trimmomatic PE -threads 5 -phred33 $RAW/${a}_1.fq.gz $RAW/${a}_2.fq.gz -baseout $TRIM/${a}_trimmed.fq ILLUMINACLIP:$ADAPT/TruSeq2OV$ LEADING:3 TRAILING:3 MINLEN:36
done
```
Details:
- -phred33: the encoding scheme used for representing base quality scores in the input sequencing data files
- -baseout: the path and name of the output files
- ILLUMINACLIP: the path to the adapters files
- LEADING: any bases at the start of a read (5') with a quality score less than 3 will be trimmed off
- TRAILING: bases at the end of a read (3') with a quality score less than 3 will be trimmed off
- MINLEN: reads with a length below 36 bases will be removed from the dataset

## Run fastQC and multiQC post-trimming using paired files to re-evaluate the content of adapters ##

# Step 3: Split of the mitogenome from the nuclear genome - BBsplit

BBsplit was designed to separate sequencing reads into different bins based on the reference databases provided by the user.

Script:
``` bash
#------Set relative path--------------------
TRIM=/home/ralvarezv/trimmomatic_01  #Directory where trimmed reads are located (input)
MITOREF=/home/ralvarezv/bbsplit_03  #Directory where the reference mitogenome is located
SPLIT=/home/ralvarezv/bbsplit_03  #Directory where splitted reads will be located (output)
#------Command------------------------------
cat list.txt | while read a
do
  bbsplit.sh ref=$MITOREF/mitogenome.fasta in1=$TRIM/${a}_trimmed_1P.fq in2=$TRIM/${a}_trimmed_2P.fq basename=${a}_mt.fq out1=$SPLIT/${a}_nomt_R1.fq out2=$SPLIT/${a}_nomt_R2.fq
done
```
Details:
- in1= path to reads (forward)
- in2= path to reads (reverse)
- basename= path and name of the file where mitogenome reads will be located (output name)
- out1= path and name of the nuclear reads (forward)
- out2= path and name of the nuclear reads (reverse)

Move the splitted mitogenomes and reference mitogenome to a new folder called "mitogen_03".

# Step 4: Map reads to the reference genome - BWA  (Burrows-Wheeler Aligner)

BWA is a software package for mapping low-divergent sequences against a large reference genome.

### Step 4a: Genome indexing
Before aligning the reads to the reference genome, it is necessary to index the reference genome. 
The reference genome used was: GCF_015237465.2

Example command:
`bwa index reference.fasta`

Script:
```
#------Set relative path--------------------
REF=/home/ralvarezv/reference  #Directory where reference genome is located (input)
#------Command------------------------------
bwa index -p $REF/cmydas_index $REF/GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna
```
Details:
- p: output database prefix
- Results should be located in the “reference” folder and they should be 5 files:
    - cmydas_index.amb
    - cmydas_index.ann
    - cmydas_index.bwt
    - cmydas_index.pac
    - cmydas_index.sa
  
These files should be moved to a new folder called "bwa_04".

### Step 4b: Genome aligment
To align our reads against the reference genome we used bwa mem.
It is the primary and most commonly used algorithm in BWA for aligning high-quality sequencing reads, particularly from next-generation sequencing (NGS) platforms like Illumina.

Here is a basic example of how to use bwa mem: `bwa mem reference.fasta reads.fq > aligned_reads.sam` 

Script:
```
#------Set relative path--------------------
REF=/home/ralvarezv/bwa_04  #Directory where reference/indexed genome is located
OUT=/home/ralvarezv/genome_aligned  #Directory where output files will be located (output)
SPLIT=/home/ralvarezv/bbsplit_03  #Directory where nuclear reads are located (input)
#------Command------------------------------
cat list.txt | while read a
do
  bwa mem -t 22 -M -R @RG\tID:${a}\tLB:CKDL230023566\tPL:ILLUMINA\tPU:A00742\tSM:${a} \
  $REF/cmydas_index $SPLIT/${a}_nomt_R1.fq $SPLIT/${a}_nomt_R2.fq > $OUT/${a}.algn.sam
done
```
Details:
- -t 22: we used 22 threads to run this process
- -M: mark shorter splits as secondary (for Picard compatibility)
- -R: complete with sample names, library IDs, platform IDs and sample IDs. They all came from the same bookstore, platform, etc. The only thing that varied was the name of the sample

# Step 5: Samtools
Samtools provides a set of utilities that allow users to manipulate and analyze data in SAM (Sequence Alignment/Map) and BAM (Binary Alignment/Map) formats.

### 5a: Convert .sam to .bam
The samtools view command converts between SAM and BAM formats. 
SAM is a human-readable text format, while BAM is a binary version that is more compact and faster to process.

Script:
```
#------Set relative path--------------------
SAM=/home/ralvarezv/genome_aligned  #Directory where the sam files are located (input)
#------Command------------------------------
cat list.txt | while read a
do
  samtools view -q 10 -f 0x2 -bSh -@ 5 $SAM/${a}.algn.sam > $BAM/${a}.algn.bam
done
```
Details:
- -q 10: skip alignments with MAPQ less than 10
- -f 0x2: FLAG, only generates alignments with all the bits configured in FLAG and checks the list of FLAGS.
          0x2: Proper pair, each segment aligns properly with its pair.
- -b: output is delivered in .sam format
- -S: automatically detect input
- -h: include the header in the output

### 5b: Check statistics
The samtools flagstat command is used to generate simple statistics from a BAM file. 
These statistics provide information about the number and types of reads in the alignment file based on their alignment flags. 
The output includes details about how many reads are properly paired, how many are singletons, and various other categories.

Script:
```
#------Set relative path--------------------
BAM=/home/ralvarezv/genome_aligned  #Directory where the bam files are located (input)
METDIR=/home/ralvarezv/genome_aligned/metrics  #Directory where the metric data will be located (output)
#------Command------------------------------
cat list.txt | while read a
do
   samtools flagstat $BAM/${a}_}.algn.bam > $METDIR/${a}_stats.txt
done
```
Details: This is an example about what we can observe in the stats.txt file for a sample:
  ```
  25781849 + 0 in total (QC-passed reads + QC-failed reads)
  22603 + 0 secondary
  0 + 0 supplementary
  0 + 0 duplicates
  25781849 + 0 mapped (100.00% : N/A)
  25759246 + 0 paired in sequencing
  12883130 + 0 read1
  12876116 + 0 read2
  25759246 + 0 properly paired (100.00% : N/A)
  25759246 + 0 with itself and mate mapped
  0 + 0 singletons (0.00% : N/A)
  0 + 0 with mate mapped to a different chr
  0 + 0 with mate mapped to a different chr (mapQ>=5)
  ```

### 5c: Sort according to genome coordinates
The samtools sort command is used to sort a BAM file according to genomic coordinates. 
Sorting is a necessary step before many downstream analyses, as it allows for efficient and ordered access to the aligned reads.

Script:
```
#------Set relative path--------------------
SORT=/home/ralvarezv/genome_aligned/sorted   #Directory where the bam sorted files will be located (output)
BAM=/home/ralvarezv/genome_aligned/metrics   #Directory where the bam files are located (input)
#------Command------------------------------
cat list.txt | while read a
do
  samtools sort -o $SORT/${a}.algn.sort.bam -@ 5 $BAM/${a}.algn.bam
done
```
Details:
- -o: leave the results in the sorted folder with this name
- -@ 5: number of cpu to use

### 5d: Index reference genome
The samtools faidx command is used to create an index for a FASTA-formatted reference genome. 
This index allows for efficient retrieval of sequences from specific genomic regions, providing quick access to the underlying nucleotide information. 
The index file created by samtools faidx has the extension ".fai."

Command:
```
samtools faidx GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna
```

# Step 6: Clean bam files

### 6a: Identify and mark duplicates - Picard
Picard MarkDuplicates is a command from the Picard Tools suite, which is a collection of command-line tools for manipulating high-throughput sequencing data. 
The MarkDuplicates tool is specifically designed to identify and mark duplicate reads in BAM files, that may have originated from the same original DNA fragment. 
These duplicates can arise during library preparation or sequencing and may impact downstream analyses, particularly in variant calling or other applications where unique reads are essential.
The metrics file (metrics.txt) provides information about the number and types of duplicates found during the process.

Script:
```
#------Set relative path--------------------
SORT=/home/ralvarezv/genome_aligned/sorted   #Directory where the bam sorted files are located (input)
DEDUP=/home/ralvarezv/genome_aligned/deduplicated   #Directory where the bam output (deduplicated) files will be located (output)
#------Command------------------------------
cat list.txt | while read a
do
  java -jar picard.jar MarkDuplicates \
  -I $SORT/${a}.algn.sort.bam -O $DEDUP/${a}.dedup.bam -METRICS_FILE $DEDUP/dedup.metrics.txt \
  -VALIDATION_STRINGENCY LENIENT -CREATE_INDEX true -CREATE_MD5_FILE true -TAGGING_POLICY All -ASSUME_SORT_ORDER coordinate
done
```
Details:
- METRICS_FILE: .txt file with the deduplication metrics
- VALIDATION_STRINGENCY=LENIENT: Picard is less strict when validating data, will be permissive in terms of input validation.
- CREATE_INDEX true: creates an index file for the deduplicated BAM file.
- CREATE_MD5_FILE true: creates an MD5 checksum file for the deduplicated BAM file.
- TAGGING_POLICY All: Specifies the tagging policy for duplicate reads. In this case, it's set to "All," which means all duplicate reads will be marked.
- ASSUME_SORT_ORDER coordinate: Assumes that the input BAM file is sorted in coordinate order. It is important to correctly specify the sorting order to ensure accurate duplicate marking.

### Clean the BAM files - Picard
The SamClean tool of Picard cleans the provided SAM/BAM files, soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads.
```
#------Set relative path--------------------
OUT=/home/ralvarezv/genome_aligned/clean_bam   #Directory where the bam cleaned files will be located (output)
DEDUP=/home/ralvarezv/genome_aligned/deduplicated   #Directory where the deduplicated bam files are located (input)
#------Command------------------------------
cat list.txt | while read a
do
  java -jar picard.jar CleanSam \
  -I $DEDUP/${a}.dedup.bam \
  -O $OUT/${a}cleaned.bam
done
```

# Step 7 



