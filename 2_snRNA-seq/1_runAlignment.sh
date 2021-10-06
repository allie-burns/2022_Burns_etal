## #############################################################################
## Date:        February 2019
## Author:      Allison M. Burns
## Filename:    /2_snRNA-seq/1_runAlignment.sh
## Project:     Epigenetic Priming - snRNA-seq analysis
## Description: This script was used to align fastq files from the snRNA-seq
##              experiment to the mouse genome.  
## #############################################################################

## Align HDACi-CFC Treatment, replicate 1
./CellRanger/cellranger-3.0.1/cellranger-cs/3.0.1/bin/count \
    --id=Allruns_HDACi_1_CR56_expect4000 \
    --transcriptome=./CellRanger/mm10-3.0.0_premrna/ \
    --fastqs=/fastdrive/CellRangerFastqOutputFD/fastq_HDACi_1/nxid7713,/fastdrive/CellRangerFastqOutputFD/fastq_HDACi_1/nxid7750,/fastdrive/CellRangerFastqOutputFD/fastq_HDACi_1/nxid7803,/fastdrive/CellRangerFastqOutputFD/fastq_HDACi_1/nxid7809 \
    --sample=HDACi_1 \
    --expect-cells=4000 \
    --chemistry=SC3Pv3 \
    --r2-length=56

## Align HDACi-CFC Treatment, replicate 2
./CellRanger/cellranger-3.0.1/cellranger-cs/3.0.1/bin/count \
    --id=Allruns_HDACi_2_CR56_expect4000 \
    --transcriptome=./CellRanger/mm10-3.0.0_premrna/ \
    --fastqs=/fastdrive/CellRangerFastqOutputFD/fastq_HDACi_2/nxid7713,/fastdrive/CellRangerFastqOutputFD/fastq_HDACi_2/nxid7750,/fastdrive/CellRangerFastqOutputFD/fastq_HDACi_2/nxid7803,/fastdrive/CellRangerFastqOutputFD/fastq_HDACi_2/nxid7809 \
    --sample=HDACi_2 \
    --expect-cells=4000 \
    --chemistry=SC3Pv3 \
    --r2-length=56

## Align Veh-CFC Treatment, replicate 1
./CellRanger/cellranger-3.0.1/cellranger-cs/3.0.1/bin/count \
    --id=Allruns_Veh_1_CR56_expect4000 \
    --transcriptome=./CellRanger/mm10-3.0.0_premrna/ \
    --fastqs=/fastdrive/CellRangerFastqOutputFD/fastq_Veh_1/nxid7713,/fastdrive/CellRangerFastqOutputFD/fastq_Veh_1/nxid7750,/fastdrive/CellRangerFastqOutputFD/fastq_Veh_1/nxid7803,/fastdrive/CellRangerFastqOutputFD/fastq_Veh_1/nxid7809 \
    --sample=Veh_1 \
    --expect-cells=4000 \
    --chemistry=SC3Pv3 \
    --r2-length=56

## Align Veh-CFC Treatment, replicate 2
./CellRanger/cellranger-3.0.1/cellranger-cs/3.0.1/bin/count \
    --id=Allruns_Veh_2_CR56_expect4000 \
    --transcriptome=./CellRanger/mm10-3.0.0_premrna/ \
    --fastqs=/fastdrive/CellRangerFastqOutputFD/fastq_Veh_2/nxid7713,/fastdrive/CellRangerFastqOutputFD/fastq_Veh_2/nxid7750,/fastdrive/CellRangerFastqOutputFD/fastq_Veh_2/nxid7803,/fastdrive/CellRangerFastqOutputFD/fastq_Veh_2/nxid7809 \
    --sample=Veh_2 \
    --expect-cells=4000 \
    --chemistry=SC3Pv3 \
    --r2-length=56
