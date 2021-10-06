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

perl ./4_rmDupByMids.pl -mode="all" \
    -in "/work/upgraeff/Allie/2_ChIP_context/3_bowtieAlign" \
    -out "/work/upgraeff/Allie/2_ChIP_context/4_bam_deDup" \
    -samtools "/ssoft/spack/humagne/v1/opt/spack/linux-rhel7-x86_S6g1_Mellanox/gcc-7.4.0/samtools-1.9-armssgpe75ofxslcuzcaobarwqgtt7o2/bin/"
