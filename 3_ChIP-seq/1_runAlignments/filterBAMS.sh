## Filter BAMS
samtools view -q 40 -b ./nxid14957_Veh1_input_S1_sorted.rmdup.bam > ../ChIP_rm_multi/nxid14957_Veh1_input_S1_sorted.rmdup.bam

samtools view -q 40 -b ./ChIP_rmdup_bams/nxid14958_Veh1_H3K27ac_S2_sorted.rmdup.bam > ./ChIP_rm_multi/nxid14958_Veh1_H3K27ac_S2_sorted.rmdup.bam 

samtools view -q 40 -b ./ChIP_rmdup_bams/nxid14959_Veh2_input_S3_sorted.rmdup.bam > ./ChIP_rm_multi/nxid14959_Veh2_input_S3_sorted.rmdup.bam

samtools view -q 40 -b ./ChIP_rmdup_bams/nxid14960_Veh2_H3K27ac_S4_sorted.rmdup.bam > ./ChIP_rm_multi/nxid14960_Veh2_H3K27ac_S4_sorted.rmdup.bam 

samtools view -q 40 -b ./ChIP_rmdup_bams/nxid14961_Veh3_input_S5_sorted.rmdup.bam > ./ChIP_rm_multi/nxid14961_Veh3_input_S5_sorted.rmdup.bam

samtools view -q 40 -b ./ChIP_rmdup_bams/nxid14962_Veh3_H3K27ac_S6_sorted.rmdup.bam > ./ChIP_rm_multi/nxid14962_Veh3_H3K27ac_S6_sorted.rmdup.bam

samtools view -q 40 -b ./ChIP_rmdup_bams/nxid14963_HDACi1_input_S7_sorted.rmdup.bam > ./ChIP_rm_multi/nxid14963_HDACi1_input_S7_sorted.rmdup.bam

samtools view -q 40 -b ./ChIP_rmdup_bams/nxid14964_HDACi1_H3K27ac_S8_sorted.rmdup.bam > ./ChIP_rm_multi/nxid14964_HDACi1_H3K27ac_S8_sorted.rmdup.bam 

samtools view -q 40 -b ./ChIP_rmdup_bams/nxid14965_HDACi2_input_S9_sorted.rmdup.bam > ./ChIP_rm_multi/nxid14965_HDACi2_input_S9_sorted.rmdup.bam

samtools view -q 40 -b ./ChIP_rmdup_bams/nxid14966_HDACi2_H3K27ac_S10_sorted.rmdup.bam > ./ChIP_rm_multi/nxid14966_HDACi2_H3K27ac_S10_sorted.rmdup.bam

samtools view -q 40 -b ./ChIP_rmdup_bams/nxid14967_HDACi3_input_S11_sorted.rmdup.bam > ./ChIP_rm_multi/nxid14967_HDACi3_input_S11_sorted.rmdup.bam

samtools view -q 40 -b ./ChIP_rmdup_bams/nxid14968_HDACi3_H3K27ac_S12_sorted.rmdup.bam > ./ChIP_rm_multi/nxid14968_HDACi3_H3K27ac_S12_sorted.rmdup.bam

## Index BAMS
samtools index ./ChIP_rmdup_bams/*.bam
