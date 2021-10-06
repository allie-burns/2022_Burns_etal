#!/bin/bash -l

## #############################################################################
## Date:        December 2020
## Author:      Allison M. Burns
## Filename:    /3_ChIP-seq/2_trimAdapt.sh
## Project:     Epigenetic Priming - ChIP-seq analysis
## Description: Trim adapters (from library preparation) and output updated
##              fastq files. This file is made to be run on a Slurm cluster.
## #############################################################################

#SBATCH -J cutAdapt
#SBATCH --nodes 1-1
#SBATCH --exclusive
#SBATCH --mem 60000
#SBATCH --time 48:00:00
  
module purge 
module load trimmomatic

## Get list of files
path="."
samples=$(ls $path/1_MIDappend_fastq/*.fastq.gz | sed 's#.*/##' | rev | cut -c10- | rev)

for samp in $samples
do
    ## Set up file system for STAR runs
    echo Trimming: $samp

    ## Get fastq files
    fastq=$(ls $path/1_MIDappend_fastq/*.fastq.gz | grep -e $samp ) 

    ## Get adapter fasta file
    adapter_seq=$(ls $path/trimAdapters/*.fasta | grep -e $samp ) 

    ## Run trimmomatic
    trimmomatic SE -phred33 \
        -threads 8 \
	$fastq \
        $path/2_trimmed_fastq/${samp}_trimmed.fastq.gz \
        -trimlog $path/2_trimmed_fastq/trim_logs/${samp}_trimlog \
        -summary  $path/2_trimmed_fastq/summary_logs/${samp}_summary \
        ILLUMINACLIP:$adapter_seq:0:6:6 \
	SLIDINGWINDOW:10:20 \
	MINLEN:36

done
