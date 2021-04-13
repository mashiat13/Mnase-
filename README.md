# Mnase-seq
We are interested in looking at histone retention in the sperm of WT and mutant mice with lysine 49 residue altered to alanine. In order to do that we performed Mnase Seq, where we added ecoli DNA in 1:1000 ratio to normalize. 

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
Then I moved on to running the alignment through Fseq2 pipeline. I was having format recognition problems with my Bam files from the alignment. So I converted them to bed files using bedtools 'bamtobed' and then ran the Fseq2 pipeline. For the control input DNA, I am using a ChIP seq sequence from WT sperm input DNA. 

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
I used a python script to normalise the np_array .h5 signal file. So far the scaling ratios I have used are 0.001, 0.004 and 0.009. (I think this definitely has to be changed). The python script is as follows:

```
import numpy as np
import h5py
import pyBigWig
from scipy.ndimage import filters
from scipy.signal import fftconvolve

def gaussian_kernel_fft(sig, sigma):
    lw = int(4 * sigma + 0.5)

    return fftconvolve(sig, filters._gaussian_kernel1d(sigma, 0, lw), mode='same')
    
def read_chrom_size_file(file_name):
    with open(file_name, 'r') as file:
        chr_size_dic = {}
        for line in file:
            line = line.strip().split()
            chr_size_dic[line[0]] = int(line[1])

    return chr_size_dic
    
chrom_ls = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10','chr11', 'chr12','chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19','chrX', 'chrY'] # chrs to deal with
#chrom_ls = ['chr1', 'chr2'] # chrs to deal with
ecoli_human_ratio = 0.001 # change the scale #0.001 number of mapped reads to ecoli/number of mapped reads to mm10

sigma_value = int(1000/6) ## INCREASE SIGMA FOR SMOOTHER SIGNAL ##

# read in sig
with h5py.File('Sample_3067-MR-1_fseq_sig.h5', mode='r') as sig_file: #change name accordingly
    sig_dict = {}
    for chrom in chrom_ls:
        signal = sig_file[chrom][:]
        signal = signal * ecoli_human_ratio ## SCALE HERE ##
        signal = np.round(gaussian_kernel_fft(signal, sigma=sigma_value), 2).astype(np.float16) ## SMOOTH HERE ##
        sig_dict[chrom] = signal
    
    # output is a signal dictionary which is indexed by chroms
    # notice array starting locations are not included, use sig_file.attrs[chrom] to check
    # if wants to output bigwig, run the following:
    
    # read in chrom_size which is needed for output bigwig
    chr_size_dic = read_chrom_size_file('/scratch/junzli_root/junzli/mrabbani/3067-MR/fastqs_3067-MR/mm10.chrom.sizes')

    # output sig in bigwig format
    with pyBigWig.open('scaled_all_chromosomes_MR1_001.bw', 'w') as output_bw: #change output name accordingly
        output_bw.addHeader([(chrom, chr_size_dic[chrom]) for chrom in chrom_ls])
        for chrom in chrom_ls:
            output_bw.addEntries(chrom, int(sig_file.attrs[chrom]), values=sig_dict[chrom], span=1, step=1)
```

I looked at the generation of the raw and the scaled bigwig files using deeptools:

![image](https://user-images.githubusercontent.com/54853508/114576750-6c647280-9c49-11eb-86d7-c102a17c4571.png)

Correlation between samples were increased from ~60% to ~80% after scaling. 

The script I used to generate the plots is as follows:
```
# The application(s) to execute along with its input arguments and options
module load Bioinformatics
module load python/3.7.4

#generate big wig files using deeptools

/home/mrabbani/.local/bin/multiBigwigSummary bins -b Sample_3067-MR-1/fseq_output/scaled_all_chromosomes_MR1_001.bw Sample_3067-MR-3/fseq_output/scaled_all_chromosomes_MR3_001.bw -o mnase_wt_mt_multibw_001.npz

/home/mrabbani/.local/bin/plotCorrelation --corData mnase_wt_mt_multibw_001.npz --corMethod pearson --whatToPlot heatmap -o mnase_wt_mt_multibw_corr_pearson_heatmap_001.pdf

/home/mrabbani/.local/bin/plotPCA -in mnase_wt_mt_multibw_001.npz -o mnase_wt_mt_multibw_PCA_readCounts_001.png -T "PCA of Mnase of WT vs Mutant"
```

I also looked at the heatmap of the reads around the Transcription start sites, using deeptools. For this I used the coding exon track from UCSC.
MR1_scaled_001_codingExon_TSS_Heatmap.png![image](https://user-images.githubusercontent.com/54853508/114583148-4b068500-9c4f-11eb-8561-1b98d1d17312.png)MR3_scaled_001_codingExon_TSS_Heatmap.png![image](https://user-images.githubusercontent.com/54853508/114583163-4fcb3900-9c4f-11eb-86e6-baf65839502b.png)

Script used to generate plots are as follows:

```
# The application(s) to execute along with its input arguments and options
module load Bioinformatics
module load python/3.7.4

/home/mrabbani/.local/bin/computeMatrix reference-point --referencePoint TSS -b 10000 -a 10000 -R bed_ucsc/mm10_ucsc_coding_exons.bed -S Sample_3067-MR-1/fseq_output/scaled_all_chromosomes_MR1_001.bw --skipZeros -o deeptools_output/MR1_scaled_001_codingExon_TSS.gz --outFileSortedRegions regions_MR1_scaled_001_codingExon_TSS.bed -p 8

/home/mrabbani/.local/bin/computeMatrix reference-point --referencePoint TSS -b 10000 -a 10000 -R bed_ucsc/mm10_ucsc_coding_exons.bed -S Sample_3067-MR-3/fseq_output/scaled_all_chromosomes_MR3_001.bw --skipZeros -o deeptools_output/MR3_scaled_001_codingExon_TSS.gz --outFileSortedRegions regions_MR3_scaled_001_codingExon_TSS.bed -p 8

/home/mrabbani/.local/bin/plotHeatmap -m MR1_scaled_001_codingExon_TSS.gz -out MR1_scaled_001_codingExon_TSS_Heatmap.png 

/home/mrabbani/.local/bin/plotHeatmap -m MR3_scaled_001_codingExon_TSS.gz -out MR3_scaled_001_codingExon_TSS_Heatmap.png 

```
In addtion to using Fseq2, I also used MACS2 to call peaks:


I used the ChipSeeker package on R to look at distribution of peaks across different genomic regions:


