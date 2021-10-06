## ##############################################################################
## Date:        March 2021
## Author:      Allison M. Burns
## Filename:    /1_RNA-seq/3_DEanalysis.R
## Project:     Epigenetic Priming - RNA-seq analysis
## Description: Run the differential expression analysis (for the chosen brain
##              region).  Inputs include sample information table (sampInfo.rds)
##              and Gene Counts for all samples (GeneCounts.rds).  This script
##              uses DESeq2 to run the Wald test for the comparisons defined
##              (all groups vs Veh-Context (VC) and the HDACi comparisons between
##              CFC (HS) and Context (HC). It then saves the output of the DE
##              analysis in R readable format and as an xlsx file.
## ##############################################################################

library(DESeq2)
library(gprofiler2)

################################################################################
## Setup
################################################################################
## Set brain region of interest
BR <- "hpp" ## "str"

## load sample information
sampInfo <- readRDS("./files/sampInfo.rds")

## Load gene counts
counts <- readRDS("./files/GeneCounts.rds")
counts <- do.call(cbind,counts)

## Make sure counts and sampInfo are in the same order
identical(colnames(counts),sampInfo$sample)

## Setup Summarized Count data 
counts <- DESeqDataSetFromMatrix(counts,
                                 colData = sampInfo,
                                 design = ~group)

## Filter counts by brain region of interest
counts <- counts[,colData(counts)$brain == BR ]

## Filter any genes that don't have atleast one count in 4 samples
keep <- rowSums(counts(counts) >1 ) >= 4 ## I chose 4 because # of replicates
counts <- counts[keep,] ## 

################################################################################
## Differential Expression Analysis
################################################################################
countsMF <- counts  ## copy counts for possible re-run later

## Set up DE comparisons
countsMF <- DESeq(countsMF,
                  test = "Wald",
                  fitType = "local")

comparisons <- list(c("VC","VS"),
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
        file = paste("./files/DEseq2_", BR,".rds", sep = ""))


writexl::write_xlsx(lapply(res,data.frame),
                    path = paste("./files/DEseq2_", BR,".xlsx",
                                 sep = "")
                    )
