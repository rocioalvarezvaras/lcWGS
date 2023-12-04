# lcWGS GREEN SEA TURTLE

---

# UPSTREAM ANALYSIS

# Step 1a: QC raw reads-fastQC

- Input: fq.gz files
- Output: .html and .zip files per sample (forward and reverse, separately)

Execute script: ./step_1.sh

# Step 1b: QC raw reads-multiQC

- Input: .html files all samples (forward and reverse, separately)
- Output: one .html file (all samples together)

Execute multiqc: $multiqc .

# Step 2: Trimming raw reads-trimmomatic

- Input: fq.gz files
- Output: four files per sample, two unpaired (U) and two paired (P)
- Add list of adapters (Illumina)

Execute script: ./step_2.sh

## Run fastQC and multiQC post-trimming using paired files to re-evaluate the content of adapters
