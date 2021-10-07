## #############################################################################
## Date:        December 2020
## Author:      Allison M. Burns
## Filename:    /3_ChIP-seq/7_ChromHMM_Homer.sh
## Project:     Epigenetic Priming - ChIP-seq analysis
## Description: Run chromHMM on bam files for multiple histone post-translational
##              modifications (in our case, we used data coming from Halder et
##              al. 2016 collected 1 hour after CFC).  First, we need to binerze
##              the bam files for the mm10 genome then we can run the chromatin
##              model (allowing for 8 chromatin states).  Then run Homer
##              annotatePeaks.pl in order to assign peaks from DiffBind to gene
##              regions based on proximity. 
## #############################################################################

################################################################################
## Command for ChromHMM
################################################################################
## Binerize BAMS
java -mx1600M -jar ./ChromHMM/ChromHMM.jar BinarizeBam -gzip ./ChromHMM/CHROMSIZES/mm10.txt ./Halder_Downloads/0_alignment/filtered_bams ./cellmarkerfiletable.txt ./BinBamDir

## Learn chromatin Model
java -mx1600M -jar ./ChromHMM/ChromHMM.jar LearnModel -p 2 ./BinBamDir ./LearnModelDir_8 8 mm10

################################################################################
## Command for Homer
################################################################################
annotatePeaks.pl homer_input_peaks.bed mm10 homer_out_peaks.txt
