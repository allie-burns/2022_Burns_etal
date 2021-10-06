#!/bin/bash

## Get list of file
path="../bams/5_bam_mapq"
samples=$(cat ../sampInfo.txt | cut -f 1,3 | awk '{print $1"_"$3""$2}' | uniq)

for samp in $samples
do
    ## Set up file system for STAR runs
    echo PeakCalling: $samp

    H3K27ac=$(ls $path/*.bam | sed 's#.*/##' | grep $samp | grep "H3K27ac" )
    input=$(ls $path/*.bam | sed 's#.*/##' | grep $samp | grep "input")

    macs2 callpeak \
	 -t $path/$H3K27ac \
	 -c $path/$input \
	 -f BAM \
	 --outdir ../peakcalls/ \
	 -n $samp \
	 --broad

done


