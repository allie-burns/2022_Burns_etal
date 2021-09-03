#!/bin/bash

#SBATCH --job-name ChIP-Align
#SBATCH --nodes 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 50G
#SBATCH --time 48:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=allison.burns@epfl.ch

module load gcc
module load bowtie2
module load samtools 

##########################################################
## STEP 1: SET UP FILE SYSTEM
##########################################################
J=ChIP
## 1.1 Make output directory
if [ ! -d ./output_$J ]; then
    mkdir ./output_$J
fi
cd ./output_$J

## 1.2 Make downstream directories within output directory
fileList=(alignmentFiles)

for f in ${fileList[@]}
do
    if [ ! -d ./$f ]; then
 	mkdir ./$f
    fi
done
 
# ##########################################################
# ## STEP 2: RUN ALIGNMENTS and SORT AND INDEX BAM FILES
# ##       - Remember to scp star_idx first
# ##       - check correct genomeDIR
# ##########################################################
samples=$(ls /work/upgraeff/Allie/3_Allie_ChIP/trimmed_fastq/*.fastq.gz | sed 's#.*/##' )


for samp in $samples
do
    ## Set up file system for STAR runs
    samp=${samp::-9}
    echo Aligning: $samp 
    fbname=$(basename $samp .fastq.gz)
    ##    mkdir ../output_$J/alignmentFiles/$fbname

    files=$(ls /work/upgraeff/Allie/3_Allie_ChIP/trimmed_fastq/*.fastq.gz | grep $samp) 

    ## Run Bowtie2
    bowtie2 -p 10 -q \
	-x /work/upgraeff/genomes/mm10_bowtie2_build/GrCm38 \
	-U $files \
	-S ./${samp}_unsorted.sam
    
    ## SAM2BAM
    samtools view -h -S -b \
	-o ./${samp}_unsorted.bam \
	./${samp}_unsorted.sam

    ## sort bam
    samtools sort ${samp}_unsorted.bam \
	-o ${samp}_sorted.bam
done

