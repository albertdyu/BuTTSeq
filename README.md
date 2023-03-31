# Description

A Snakemake pipeline for the analysis of data produced by Butt-Seq library preps. It processes paired-end Butt-Seq fastq files and returns several files that may be of interest to the investigator, including:

1. Deduplicated BAM files, with and without small RNAs computationally removed
2. Single-nucleotide resolution BAM files
3. Stringtie assemblies and counts files for potential transcript discovery
4. A featurecounts table with raw counts
5. Normalized full-read and single nucleotide resolution BW files. 

This pipeline produces analyses and outputs that not every investigator will find a use for, but some may find interesting.

## Input files

This pipeline requires 3 fastq files, formatted as such:

* ___(sample name)_R1_001.fastq.gz - Read 1
* ___(sample name)_R2_001.fastq.gz - UMI Read
* ___(sample name)_R3_001.fastq.gz - Read 2

To produce 3 read files from a sequencing run, edit this line in RunInfo.xml:

\<Read Number="2" NumCycles="8" IsIndexedRead="Y" />

to 

\<Read Number="2" NumCycles="8" IsIndexedRead="N" />

and demultiplex using Bcl2fastq as usual. 

## Output Files

1. Deduplicated BAM files, with and without small RNAs computationally removed: found in dedup/ and subSno/
2. Single-nucleotide resolution BAM files: found in SNR/
3. Stringtie assemblies and counts files for potential transcript discovery: found in stringtie/
4. A featurecounts table with raw counts: found in results/counts.featureCounts
5. Normalized full-read and single nucleotide resolution BW files, separated by strand: found in bws/paired and bws/SNR, respectively

## Contents 

- 'beds/': A bed file describing the 3' end of small RNAs and one bed file describing the 3' end of exons in the Drosophila Melanogaster assembly dm6. 

- 'scripts/': A series of scripts used to process the data, including:

    - removeclipping.py: A script from NGSUtils, slightly modified to suit this pipeline (Breese et al, 2013). This script removes softclipped reads from the BAM file prior to conversion into single nucleotide reads. Without removing softclipping first, the following script will erroneously assign the 3' most end as the soft-clipped base. 

    - get_SNR_bam.py: A script by Tomás Gomes, slightly modified to suit this pipeline (Nojima et al, 2015). Converts a bam file to only contain the first base of Read 2. 

    - deseq2_normalization.R: An Rscript to load the featurecounts table into DESeq2 and output normalization factors, which are using as scaling factors by bamCoverage when converting to BigWig files for visualization and/or metagene plotting.  

# Installation and usage

## Installation

If you are new to conda and/or snakemake, you'll need to install the appropriate version of miniconda3 for your operating system. 

1. Clone or download this github repo into a directory containing your RNA-Seq files. 
2. Inside your terminal, create a conda environment with the necessary dependencies by running the following command:


    conda env create --name butt --file envs/environment.yaml


3. Activate your new conda environment with the following command:


    conda activate butt


4. Set up your config.yaml, as described below.

5. Run the Snakemake pipeline with the following command (Alter the number of cores as desired):


    snakemake --cores 8 --configfile config.yaml


## Configuration

'config.yaml' must be edited to the following parameters:

samples:

    (SampleName_1): Path/To/Sample1
    
    (SampleName_2): Path/To/Sample2
    
    etc. etc.
    
annotation:

    Path/To/GTF/Annotation
    
index:

    Path/To/STAR/Index
    
smallRNAs:

    Path/To/Undesirable/SmallRNAs
    
exonends:

    Path/To/ExonEnds
    

Small RNAs consist of chromatin-associated RNAs that are not thought to be products of active transcription and are typically not an analyte of interest, so they are computationally removed. 

3' End ligation techniques, including Butt-Seq, often capture splicing intermediates which cannot be distinguished from true polymerase pause sites, so any reads mapping precisely to the 3' end of exons are removed. 
