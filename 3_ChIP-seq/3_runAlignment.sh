#!/bin/bash

## #############################################################################
## Date:        December 2020
## Author:      Allison M. Burns
## Filename:    /3_ChIP-seq/3_runAlignment.sh
## Project:     Epigenetic Priming - ChIP-seq analysis
## Description: Align trimmed fastq files to mm10 genome with bowtie2 then
##              convert sam to bam files and sort bam files.
## #############################################################################

#SBATCH --job-name ChIP-Align
#SBATCH --nodes 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 50G
#SBATCH --time 72:00:00

module load gcc
module load bowtie2
module load samtools 


## Get list of files
path="."
samples=$(ls $path/2_trimmed_fastq/*.fastq.gz | sed 's#.*/##' | rev | cut -c18- | rev)

for samp in $samples
do
    ## Set up file system for STAR runs
    echo Aligning: $samp 
    files=$(ls $path/2_trimmed_fastq/*.fastq.gz | grep $samp)

    ## Run Bowtie2
    bowtie2 -p 24 -q \
	-x /work/upgraeff/genomes/mm10_bowtie2_build/GrCm38 \
	-U $files \
	-S $path/3_bowtieAlign/${samp}_unsorted.sam
    
    ## SAM2BAM
    samtools view -h -S -b \
	-o $path/3_bowtieAlign/${samp}_unsorted.bam \
	$path/3_bowtieAlign/${samp}_unsorted.sam

    ## sort bam
    samtools sort $path/3_bowtieAlign/${samp}_unsorted.bam \
	-o $path/3_bowtieAlign/${samp}_sorted.bam

    rm $path/3_bowtieAlign/${samp}_unsorted.sam
    rm $path/3_bowtieAlign/${samp}_unsorted.bam 
done

