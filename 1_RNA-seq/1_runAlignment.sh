#!/bin/bash

## #############################################################################
## Date:        February 2021
## Author:      Allison M. Burns
## Filename:    /1_RNA-seq/1_runAlignment.sh
## Project:     Epigenetic Priming - RNA-seq analysis
## Description: This script was used to align fastq files from the RNA-seq
##              experiment to the mouse genome.  This code creates the file
##              system, runs the alignment and sorts and indexes the output
##              bam files. 
## #############################################################################

## #############################################################################
## Setup File System
## #############################################################################
## Output directory
if [ ! -d ./alignment ]; then
    mkdir ./alignment
fi
cd ./alignment

fileList=(STARalign STARbam)
for f in ${fileList[@]}
do
    if [ ! -d ./$f ]; then
 	mkdir ./$f
    fi
done

## #############################################################################
## Run Alignment
## #############################################################################
samples=$(cat ../sampInfo.txt | cut -f 1)

for samp in $samples ;
do
    echo Aligning: $samp 

    ## Set up file system for STAR runs
    fbname=$(basename $samp .fastq.gz)
    mkdir ../alignment/STARalign/$fbname    

    ## get files for each alignment
    files=$(ls ../fastq/*.fastq.gz | grep $samp) 

    ##Run STAR
    /work/upgraeff/Software/STAR/bin/Linux_x86_64/STAR \
      	--runThreadN 8 \
      	--genomeDir /work/upgraeff/genomes/mm10_staridx_ensemble93/ \
      	--readFilesIn $files \
    	--outFilterType BySJout \
      	--outSAMtype BAM SortedByCoordinate \
      	--alignSJDBoverhangMin 5 \
      	--alignSJoverhangMin 55 \
      	--readFilesCommand gzip -cd \
      	--outFileNamePrefix ./STARalign/$fbname/${fbname}_ &
wait

done


## #############################################################################
## Sort and index bam files
## #############################################################################
for samp in $samples ;
do
    echo Processing: $fbname

    ## Set up file system for STAR runs
    fbname=$(basename $samp .fastq.gz)

    ## Copy .bam files to STARbam
    mv ./STARalign/$fbname/${fbname}_Aligned.sortedByCoord.out.bam \
       ./STARbam/

    ## Sort and index bam files
    samtools sort -l 9 -o ./STARbam/${fbname}_sorted.bam \
	./STARbam/${fbname}_Aligned.sortedByCoord.out.bam 
    samtools index ./STARbam/${fbname}_sorted.bam

done

