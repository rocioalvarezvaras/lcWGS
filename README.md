# lcWGS GREEN SEA TURTLE

---

# UPSTREAM ANALYSIS

# Step 1a: QC raw reads-fastQC

- Input: fq.gz files
- Output: .html and .zip files per sample (forward and reverse, separately)

_Execute script: step_1.sh_

# Step 1b: QC raw reads-multiQC

- Input: .html files all samples (forward and reverse, separately)
- Output: one .html file (all samples together)

_Execute multiqc:  multiqc ._

# Step 2: Trimming raw reads-trimmomatic

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
