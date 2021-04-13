# Mnase-seq
We are interested in looking at histone retention in the sperm of WT and mutant mice with lysine 49 residue altered to alanine. In order to do that we perofmed Mnase Seq, where we added ecoli DNA in 1:1000 ratio to normalize. 

1. Fast QC of sequencing data 
2. Cutadapt to trip over-represented sequences

3. Alignment using bowtie2 

Bowtie commands:
Command for mm10 as reference 
```
# The application(s) to execute along with its input arguments and options
module load Bioinformatics
module load bowtie2/2.4.1

#this is working for us 
bowtie2 -x mm10/mm10 -1 Sample_3067-MR-1/cutadapt_3067-MR-1_CAGATCAT-AGATCTCG_S106_R1_001.fastq.gz -2 Sample_3067-MR-1/cutadapt_3067-MR-1_CAGATCAT-AGATCTCG_S106_R2_001.fastq.gz -S Sample_3067-MR-1/cutadapt_3067-MR-1_CAGATCAT-AGATCTCG_S106.sam -p 4

```
Command for ecoli genome as reference 
```
# The application(s) to execute along with its input arguments and options
module load Bioinformatics
module load bowtie2/2.4.1

#this is working for us 
bowtie2 -x ecoli/e_coli -1 Sample_3067-MR-1/cutadapt_3067-MR-1_CAGATCAT-AGATCTCG_S106_R1_001.fastq.gz -2 Sample_3067-MR-1/cutadapt_3067-MR-1_CAGATCAT-AGATCTCG_S106_R2_001.fastq.gz -S Sample_3067-MR-1/cutadapt_3067-MR-1_CAGATCAT-AGATCTCG_S106_ECOLI.sam
```
