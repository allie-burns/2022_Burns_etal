library(DiffBind)
library(rtracklayer)
library(ggplot2)

################################################################################
## Setup and read counts
################################################################################
## Load in csv file of sample information
db_obj <- dba(sampleSheet = "./data/sampleData.csv")

nPeaks <- data.frame(sample = db_obj$samples$SampleID,
                     condition = db_obj$samples$Condition,
                     treatment = db_obj$samples$Treatment,
                     group = paste(db_obj$samples$Treatment,
                                   db_obj$samples$Condition,
                                   sep = "_"),
                     nPeaks = colSums(db_obj$called))

## Calculate peak and read counts  
pdf("./figures/1_heatmap_allpeak.pdf",height = 8, width = 8)
plot(db_obj,
     main = "Occupancy (peak score) Correlation")
## Count reads - affinity analysis information
db_obj <- dba.count(db_obj)
plot(db_obj,
     main = "Affinity (read count) Correlation")
dev.off()

## Get normalized counts
db_norm <- dba.contrast(db_obj)
db_norm <- dba.analyze(db_norm, method=DBA_DESEQ2, bParallel = F, bFullLibrarySize = F)
norm_counts <- dba.report(db_norm,
                          th = 1,
                          bCounts = T,
                          DataType=DBA_DATA_FRAME,
                            bNormalized = T)

## Set contrast
db_obj <- dba.contrast(db_obj) ## 8 contrasts

## Change order of 6th contrast
db_obj$contrasts[6][[1]]$name1 <- "Context:HDACi"
db_obj$contrasts[6][[1]]$name2 <- "CFC:Veh"

g1 <- db_obj$contrasts[6][[1]]$group1
g2 <- db_obj$contrasts[6][[1]]$group2
db_obj$contrasts[6][[1]]$group1 <- g2
db_obj$contrasts[6][[1]]$group2 <- g1

################################################################################
## Differential Expression analysis - Compare EdgeR and DEseq
################################################################################
## Perform the differential analysis
db_obj <- dba.analyze(db_obj,
                      method=DBA_DESEQ2,
                      bParallel = FALSE)

myCon <- list(VS_VC = 7,
              HC_VC = 8,
              HS_VC = 5,
              HS_HC = 4,
              HC_VS = 6,
              HS_VS = 3,
              CFC_Con = 1,
              HDACi_Veh = 2)


## Make heatmaps of contrasts
lapply(names(myCon), function(name) {
    con = myCon[[name]]
    dba.plotHeatmap(db_obj, contrast = con,
                    main = paste("DE peaks - ", name))
})

## Draw volcano plots
lapply(names(myCon), function(name) {
    con = myCon[[name]]
    dba.plotVolcano(db_obj, contrast = con)
})

## Plot Overlaps
dba.plotVenn(db_obj, contrast = c(7,8,5))
dba.plotVenn(db_obj, contrast = c(1,2))
dba.plotVenn(db_obj, contrast = c(3,4,6))

## PCA plots
lapply(names(myCon), function(name) {
    con = myCon[[name]]
    dba.plotPCA(db_obj, contrast = con, label = DBA_REPLICATE)
})

## BoxPlots
lapply(names(myCon), function(name) {
    con = myCon[[name]]
    dba.plotBox(db_obj, contrast = con)
})


## Get DE reports
getReport <- function(name) {
    con = myCon[[name]]
    dba.report(db_obj, contrast = con, th = 1,method = DBA_DESEQ2)
}

report <- lapply(names(myCon),getReport)
names(report)  <-  names(myCon)

saveRDS(norm_counts, file = "./data/ChIP_DEseq2_norm_counts.rds")
saveRDS(report, file = "./data/ChIP_DEseq2_deReport.rds")
saveRDS(db_obj, file = "./data/HDACi_PeakObject.rds")
