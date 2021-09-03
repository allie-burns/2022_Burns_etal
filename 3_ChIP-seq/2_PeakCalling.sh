## Veh_1
macs2 callpeak \
-t ../bams/nxid14958_Veh1_H3K27ac_S2_sorted.rmdup.bam \
-c ../bams/nxid14957_Veh1_input_S1_sorted.rmdup.bam \
-f BAM \
--outdir ./ \
-n Veh1 \
--broad

## Veh_2
macs2 callpeak \
-t ../bams/nxid14960_Veh2_H3K27ac_S4_sorted.rmdup.bam \
-c ../bams/nxid14959_Veh2_input_S3_sorted.rmdup.bam \
-f BAM \
--outdir ./ \
-n Veh2 \
--broad

 ## Veh_3
macs2 callpeak \
-t ../bams/nxid14962_Veh3_H3K27ac_S6_sorted.rmdup.bam \
-c ../bams/nxid14961_Veh3_input_S5_sorted.rmdup.bam \
-f BAM \
--outdir ./ \
-n Veh3 \
--broad

## HDACi_1
macs2 callpeak \
-t ../bams/nxid14964_HDACi1_H3K27ac_S8_sorted.rmdup.bam \
-c ../bams/nxid14963_HDACi1_input_S7_sorted.rmdup.bam \
-f BAM \
--outdir ./ \
-n HDACi1 \
--broad

## HDACi_2
macs2 callpeak \
-t ../bams/nxid14966_HDACi2_H3K27ac_S10_sorted.rmdup.bam \
-c ../bams/nxid14965_HDACi2_input_S9_sorted.rmdup.bam \
-f BAM \
--outdir ./ \
-n HDACi2 \
--broad

## HDACi_3
macs2 callpeak \
-t ../bams/nxid14968_HDACi3_H3K27ac_S12_sorted.rmdup.bam \
-c ../bams/nxid14967_HDACi3_input_S11_sorted.rmdup.bam \
-f BAM \
--outdir ./ \
-n HDACi3 \
--broad
