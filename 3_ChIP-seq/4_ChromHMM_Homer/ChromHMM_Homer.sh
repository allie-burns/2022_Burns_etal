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
