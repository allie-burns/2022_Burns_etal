#!/bin/bash

#SBATCH --job-name bam_qualFilt
#SBATCH --nodes 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 25G
#SBATCH --time 24:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=allison.burns@epfl.ch

module load gcc
##module load bowtie2
module load samtools 

## Get list of files
path="/work/upgraeff/Allie/2_ChIP_context"
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

##samtools index $path/5_bam_mapq/*.bam
