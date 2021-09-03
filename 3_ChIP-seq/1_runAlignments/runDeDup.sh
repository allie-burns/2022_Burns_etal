#!/bin/bash

#SBATCH --job-name ChIP-deDup
#SBATCH --nodes 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 50G
#SBATCH --time 48:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=allison.burns@epfl.ch

module load gcc
module load samtools 


perl ./rmDupByMids.pl -mode="all" \
    -in "./ChIP_bams" \
    -out "./ChIP_deDup" \
    -samtools "/ssoft/spack/humagne/v1/opt/spack/linux-rhel7-x86_S6g1_Mellanox/gcc-7.4.0/samtools-1.9-armssgpe75ofxslcuzcaobarwqgtt7o2/bin/"



# samples=$(ls /work/upgraeff/Allie/3_Allie_ChIP/ChIP_alignment/*sorted.bam | sed 's#.*/##' | uniq)


# for samp in $samples
# do
#     ## Set up file system for STAR runs
#     samp=${samp::-9}
#     echo DeDupping: $samp 
#     ##    fbname=$(basename $samp .fastq.gz)
#     ##    mkdir ../output_$J/alignmentFiles/$fbname

#     files=$(ls /work/upgraeff/Allie/3_Allie_ChIP/trimmed_fastq/*.fastq.gz | grep $samp) 

#     ## Run Bowtie2
#     bowtie2 -p 10 -q \
# 	-x /work/upgraeff/genomes/mm10_bowtie2_build/GrCm38 \
# 	-U $files \
# 	-S ./${samp}_unsorted.sam
    
#     ## SAM2BAM
#     samtools view -h -S -b \
# 	-o ./${samp}_unsorted.bam \
# 	./${samp}_unsorted.sam

#     ## sort bam
#     samtools sort ${samp}_unsorted.bam \
# 	-o ${samp}_sorted.bam
# done
