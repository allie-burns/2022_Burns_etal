library(DESeq2)
library(gprofiler2)

################################################################################
## Setup
################################################################################
## load sample information
sampInfo <- readRDS("./data/sampInfo.rds")

## Load gene counts
counts <- readRDS("./data/GeneCounts.rds")
counts <- do.call(cbind,counts)

## Remove Naive groups
counts <- counts[,-grep("N", colnames(counts))]
sampInfo <- sampInfo[-grep("Nai", sampInfo$treat),]
sampInfo <- droplevels(sampInfo)

## Make sure counts and sampInfo are in teh same order
identical(colnames(counts),sampInfo$sample)

## Setup Summarized Count data 
counts <- DESeqDataSetFromMatrix(counts,
                                 colData = sampInfo,
                                 design = ~group)

## Set brain region of interest
BR <- "str"
counts <- counts[,colData(counts)$brain == BR ]

## Filter any genes that don't have atleast one count in 4 samples
keep <- rowSums(counts(counts) >1 ) >= 4 ## I chose 4 because # of replicates
counts <- counts[keep,] ## 

## Save normalized counts
dds <- estimateSizeFactors(counts)
normalized_counts <- counts(dds, normalized=TRUE)
saveRDS(normalized_counts,
        paste("./data/DE_conditioning/normCounts_",BR,".rds",sep = ""))

################################################################################
## Differential Expression Analysis
################################################################################
countsMF <- counts  ## copy counts for possible re-run later

## Set up DE comparisons
countsMF <- DESeq(countsMF,
                  test = "Wald",
                  fitType = "local")

comparisons <-  list(c("VC","VS"),
                     c("VC","HC"),
                     c("VC","HS"),
                     c("HC","HS"))

## Get result table for each comparison
getComparisonResults <- function(comparisons) {
    c1 <- comparisons[1]
    c2 <- comparisons[2]
    
    res <- results(countsMF,
                   contrast= c("group", c2, c1),
                   listValues = c(1,-1))
    
    res <- na.omit(res)
    res <- res[order(res$padj),]
    
    gene <- gconvert(rownames(res),
                     organism = "mmusculus",
                     filter_na = FALSE)
    
    res$gene  <-  gene$name
    res$ensembl_id  <- rownames(res)
    res$description  <-  gene$description
    res  <- res[,c(which(colnames(res)=="ensembl_id"),
                   which(colnames(res)!="ensembl_id"))]
    
    return(res)
}

res <- lapply(comparisons,getComparisonResults)
names(res) <- unlist(lapply(comparisons, function(x) {paste(x[[2]],x[[1]], sep = "-")}))

## save results
saveRDS(res,
        file = paste("./data/DE_conditioning/DEseq2_", BR,"_VCcontrol.rds", sep = ""))


writexl::write_xlsx(lapply(res,data.frame),
                    path = paste("./data/DE_conditioning/DEseq2_", BR,"_VCcontrol.xlsx",
                                 sep = "")
                    )
