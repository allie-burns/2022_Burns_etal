#!/bin/bash

#SBATCH --job-name RNAseq_hdac
#SBATCH --nodes 1
#SBATCH --cpus-per-task 34
#SBATCH --mem 50G
#SBATCH --time 48:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=allison.burns@epfl.ch

module load gcc
module load samtools

##########################################################
## STEP 1: SET UP FILE SYSTEM
##########################################################
# 1.1 Make output directory
if [ ! -d ./alignment ]; then
    mkdir ./alignment
fi
cd ./alignment

## 1.2 Make downstream directories within output directory
fileList=(STARalign STARbam)

for f in ${fileList[@]}
do
    if [ ! -d ./$f ]; then
 	mkdir ./$f
    fi
done

##########################################################
## STEP 2: RUN ALIGNMENTS
##########################################################
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


##########################################################
## STEP 3: sort and index bamfiles
##########################################################
for samp in $samples ;
do
    echo Processing: $fbname

    ## Set up file system for STAR runs
    fbname=$(basename $samp .fastq.gz)

    ## Copy .bam files to STARbam
    # mv ./STARalign/$fbname/${fbname}_Aligned.sortedByCoord.out.bam \
    # 	./STARbam/

    ## Sort and index bam files
    samtools sort -l 9 -o ./STARbam/${fbname}_sorted.bam \
	./STARbam/${fbname}_Aligned.sortedByCoord.out.bam 
    samtools index ./STARbam/${fbname}_sorted.bam
    ##rm ./STARbam/${fbname}_Aligned.out.bam 

done

