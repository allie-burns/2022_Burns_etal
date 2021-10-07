#!/bin/bash

## #############################################################################
## Date:        December 2020
## Author:      Allison M. Burns
## Filename:    /3_ChIP-seq/4_filterBamFiles.sh
## Project:     Epigenetic Priming - ChIP-seq analysis
## Description: This script calls the perl script (provided by Active Motif)
##              that removes duplicated reads based on MID sequences.  Be sure
##              to input bam files (aligned and sorted) and point to correct
##              samtools location. Then remove alignments with low mapping
##              quality (mapq <= 40) and re-index bam files.
## #############################################################################

## Remove duplicated reads
perl ./4_rmDupByMids.pl -mode="all" \
    -in "./files/3_bowtieAlign" \
    -out "./files/4_bam_deDup" \
    -samtools "./samtools-1.9-armssgpe75ofxslcuzcaobarwqgtt7o2/bin/"


## Remove alignments with low mapping quality
path="./files"
samples=$(ls $path/4_bam_deDup/*.bam | sed 's#.*/##' | rev | cut -c18- | rev)

for samp in $samples
do
    ## Set up file system for STAR runs
    echo Aligning: $samp 
    files=$(ls $path/4_bam_deDup/*.bam | grep $samp)

    ## Filter for mapq
    samtools view -q 40 \
	-b $files \
	> $path/5_bam_mapq/${samp}_sorted.rmdup.bam

    samtools index $path/5_bam_mapq/${samp}_sorted.rmdup.bam
done
