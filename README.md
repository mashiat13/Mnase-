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
The I moved on to running the alignment through Fseq2 pipeline. I was having format recognition problems with my Bam files from the alignment. So I converted them to bed files using bedtools 'bamtobed' and then ran the Fseq2 pipeline. For the control input DNA, I am using a ChIP seq sequence from WT sperm input DNA. 

BEDtools command to covert bam to bed files. 

```
# The application(s) to execute along with its input arguments and options
module load Bioinformatics
module load python3.6-anaconda
module load bedtools2/2.29.2

bedtools bamtobed -i Sample_3067-MR-1/bamfiles/cutadapt_3067-MR-1_CAGATCAT-AGATCTCG_S106.sorted_name.bam -bedpe > Sample_3067-MR-1/bamfiles/cutadapt_3067-MR-1_CAGATCAT-AGATCTCG_S106.sorted_name.bed 

bedtools bamtobed -i Sample_3067-MR-3/bamfiles/cutadapt_3067-MR-3_GTCCGCAT-AGATCTCG_S107.sorted_name.bam -bedpe > Sample_3067-MR-3/bamfiles/cutadapt_3067-MR-3_GTCCGCAT-AGATCTCG_S107.sorted_name.bed

bedtools bamtobed -i chip_input/cutadapt_1619-SS-12_GATTCTGC-CTCTCGTC_S16.sorted_name.bam -bedpe > chip_input/cutadapt_1619-SS-12_GATTCTGC-CTCTCGTC_S16.sorted_name.bed
```
Fseq2 command: I used the same command to output files in both np_array and bigwig format. 

```
# The application(s) to execute along with its input arguments and options
module load Bioinformatics
module load python3.6-anaconda
module load bedtools2/2.29.2

#conda create --name fseq2_env python=3.6
source activate fseq2_env
conda activate fseq2_env

fseq2 callpeak Sample_3067-MR-1/bamfiles/cutadapt_3067-MR-1_CAGATCAT-AGATCTCG_S106.sorted_name.bed -control_file chip_input/cutadapt_1619-SS-12_GATTCTGC-CTCTCGTC_S16.sorted_name.bed -l 1000 -t 4.0 -p_thr 0.05 -sig_format np_array -chrom_size_file mm10.chrom.sizes -v -cpus 4 -o Sample_3067-MR-1/fseq_output -name Sample_3067-MR-1_fseq -pe

```




