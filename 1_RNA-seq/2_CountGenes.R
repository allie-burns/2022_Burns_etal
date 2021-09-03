################################################################################
## Count reads that align to ensembl 93
################################################################################
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)

##txdb <- makeTxDbFromEnsembl(organism="Mus musculus", release=93)
##saveDb(txdb,"../annotations/mmusculus_ensembl_release93.txdb")
txdb <- loadDb("../annotations/mmusculus_ensembl_release93.txdb")

## Get rid of the unmapped, haplotypic and random chromosomal assembly
isActiveSeq(txdb) <- !grepl("_|\\.",seqnames(seqinfo(txdb)))

geneModel <- exonsBy(txdb, "gene")

## Get data from bam files to compare to annotation
bamfiles <- "../1_alignment/alignment/STARbam"
bamFiles <- BamFileList(dir(bamfiles,pattern=".bam$",recursive=TRUE,full=TRUE))

## Count gene expression data in bam files
counterPerChr <- function(idx,data,gnModel=geneModel,mapq=10){
    ## Assign chromosome and file name
    seqname <- data$chr[idx]
    fl <- data$file[idx]
    
    ## Select genes on respective chromosome
    gnModel <- geneModel[seqnames(unlist(range(geneModel))) == seqname]
    
    ## Load bamfiles
    s.i <- seqinfo(BamFile(fl)) 
    seq.length <- seqlengths(s.i)[seqname]
    
    ## Read in bam file reads
    param <- ScanBamParam(flag = scanBamFlag(
                              isProperPair = TRUE,
                              isUnmappedQuery = FALSE,
                              hasUnmappedMate = FALSE,
                              isSecondaryAlignment = FALSE,
                              isNotPassingQualityControls = FALSE,
                              isDuplicate = FALSE,
                              ),
                            which=GRanges(seqname,IRanges(1, seq.length)),
                            what='mapq',
                            tag='NH')
    aln <- readGAlignmentPairs(fl, param=param, strandMode=2)
    aln <- aln[countOverlaps(aln, gnModel)==1 ] ## No aligments overlapping >1 gene

    ## Count alignments per gene
    cat(paste("Counting overlaps for",basename(fl),"in chromosome",seqname,"\n"))
    counts <- countOverlaps(gnModel, aln)  ##Counts anything that overlaps annotated exons
    names(counts) <- names(gnModel)
    counts
}

toProcess <-
    do.call(rbind,lapply(bamFiles,function(file){ data.frame(chr=seqnames(seqinfo(file)),
                                                             file=path(file),
                                                             stringsAsFactors=FALSE)
    }))

toProcess <- toProcess[!grepl("_|\\.",toProcess$chr),]

## Reading the bam files, in parallels one chromosome at a time
counts.raw <- lapply(seq(to=nrow(toProcess)),
                     counterPerChr,
                     toProcess)

counts <- lapply(split(counts.raw,f=toProcess$file),
                 function(cnts) do.call(c,cnts))

### Just to be Kosher, reoder the genes in the result tables same as in gnModel
counts <- lapply(counts,function(x) x[names(geneModel)])
names(counts) <- sub(".bam",'',basename(names(counts)))
names(counts) <- substr(names(counts),1,nchar(names(counts))-7)

##save(counts, file=file.path(paste(dataFolder,"GeneCounts.Rda",sep="_")))
saveRDS(counts,file = "./data/GeneCounts.rds")

