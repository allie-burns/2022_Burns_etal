## #############################################################################
## Date:        December 2020
## Author:      Allison M. Burns
## Filename:    /3_ChIP-seq/1_appendMIDs.R
## Project:     Epigenetic Priming - ChIP-seq analysis
## Description: Append MID information from Read_1 to Read_2 for downstream
##              processing (and removal of duplicates).  Also make adapter
##              sequence files here for next step, adapter trimming
##              (2_trimAdapt.sh)
## #############################################################################

library(ShortRead)
library(Biostrings)

################################################################################
## Setup environment
################################################################################
## Call get sample infor for alignment
fastqdir <- "./raw_fastq"
files <- list.files(fastqdir, pattern = "fastq.gz", full.names = TRUE)

sample <- unique(substr(basename(files), 
                        1, 
                        nchar(basename(files)) - 19))

sample <- sub("_$", "", sample)

################################################################################
## Append MIDs from Read1 to Read2 
################################################################################
appendFastQs <- function(sample,files) {
    ## Get sample and fiels
    samp <- sample
    fls <- files[grep(samp, files)]

    ## Don't load full thing in at once
    r1 <- FastqStreamer(fls[grep("R1", fls)], n = 100000, readerBlockSize = 100000)
    r2 <- FastqStreamer(fls[grep("R2", fls)], n = 100000, readerBlockSize = 100000)
    on.exit(close(r1))
    on.exit(close(r2))
    
    ## Run appending
    repeat {
      ## input chunk
      r1.s <- yield(r1)
      ##cat(paste(length(r1.s),"\n"))
      r2.s <- yield(r2)
      if (length(r1.s) == 0)
        break

      ## Get strings for merging
      r1_id <- as.character(id(r1.s))
      r2_seq <- as.character(sread(r2.s))    

      r1_new <- BStringSet(
        unlist(
          lapply(seq(1:length(r1_id)), function(i){
            r1_split <- strsplit(r1_id[i], " ")
            paste(r1_split[[1]][1],"_", r2_seq[i], " ",r1_split[[1]][2],sep = "")
          })))

      ## Build new reads
      r3 <- ShortReadQ(sread(r1.s), quality(r1.s), r1_new)
      ## write the new fastq files
      writeFastq(r3,
                 paste("./fastq_append_MID/",sample,".fastq.gz", sep = ""),
                 "a")
    }
}

lapply(sample, appendFastQs,files)

################################################################################
## Make adapter files for next step (trimming adapters)
################################################################################
adapters <- read.delim("./trimAdapters/adapter_list.txt")

lapply(seq(1:nrow(adapts)), function(i) {
    sample <- adapters[i,1]
    adapt <- adapters[i,2]
    adapt <- DNAStringSet(x=adapt)
    names(adapt)  <-  sample
    writeXStringSet(adapt,
                    paste("./trimAdapters/", sample, ".fasta",sep = ""),
                    format="fasta")
})
