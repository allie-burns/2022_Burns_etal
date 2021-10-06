## #############################################################################
## Date:        April 2021
## Author:      Allison M. Burns
## Filename:    /2_snRNA-seq/2_SeuratAnalysis.R
## Project:     Epigenetic Priming - snRNA-seq analysis
## Description: This script loads and formats the count matrix from CellRanger
##              (alternative processed count matrix can be loaded at line 49),
##              makes the Seurat object, removes doublets, filters nuclei and
##              normalizes nuclear counts then clusters nuclei using TSNE and
##              UMAP clustering.
##              Then nuclear clusters are assigned to cell types based on cell-
##              type marker expression.  
## #############################################################################

library(Seurat)
library(dplyr)
library(DoubletFinder)
library(sctransform)
library(ggplot2)
library(cowplot)

## Set Plotting Colors
cols <- c("#8B829C", "#5A486F","#99CCBA", "#689D6C") 

################################################################################
## Load count matrices for all files
################################################################################
## Loading in the feature_bc_matrix from CellRanger Directory
## Files aren't included in the files dir (alternative count upload on line 49)
sDir <- list.dirs("../0_CellRangerOut/count_allruns",recursive=TRUE)
myDir <- sDir[grep("filtered_feature_bc_matrix",sDir)]
names(myDir) <- lapply(strsplit(myDir,"_"),function(x) { paste(x[4:5], collapse="_") } )
names(myDir) <- gsub("_",".",names(myDir))
myDir <- myDir[c(3,4,1,2)]

getFiles <- function(i, myDir) {
    ## Get loading information
    matrix_dir = paste(myDir[[i]],"/",sep="")
    list.files(matrix_dir)
    ##Load 10x
    data <- Read10X(data.dir = matrix_dir)
    ## rename nuclei for matrix merging 
    colnames(data) <- paste(colnames(data),names(myDir[i]),sep="_")
    data
}

count.data <- lapply(seq(1:length(myDir)),getFiles,myDir)
names(count.data) <- names(myDir)

## Merge count dataframes
counts <- do.call(cbind,count.data)
counts <- counts[!duplicated(rownames(counts)),]

## Alternatively, download the pre-processed count file
counts <- read.delim("./files/snRNA_raw_counts_matrix.txt")

################################################################################
## Create and process seurat object
################################################################################
## Set filter thresholds
min.cells <- 10 
min.genes <- 300
max.genes <- 6000 
min.UMI  <-  700
max.UMI <- 30000
max.mito <- 3  

## Setup Seurat Object
seu <- CreateSeuratObject(counts = counts, min.cells = min.cells)
seu[["percent.mt"]] <- PercentageFeatureSet(object= seu,pattern = "^mt-") ## mito counts
seu[["orig.ident"]] <- factor(sub(".*_", "", colnames(seu)), ## Rename to drug and rep #
                              levels = c("Veh.1","Veh.2","HDACi.1","HDACi.2"))
seu[["drug"]] <- factor(gsub(".[1-9]","",seu@meta.data$orig.ident),
                        levels = c("Veh","HDACi"))
seu[["day"]] <- factor(gsub("[a-zA-Z].","",seu@meta.data$orig.ident))

p1 <- VlnPlot(object = seu,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by= "orig.ident",
              pt.size = 0.25,
              cols=cols,
              ncol=3)

################################################################################
## Find and Remove Doublets
################################################################################
## Normalize Data for finding Doublets
seuDoub = seu %>%
    ## Split the object into a list for input the Seurat integration
    SplitObject(split.by = 'orig.ident') %>%
    ## Normalize the data using regularized negative binomial models
    purrr::map(~ SCTransform(.))

## Find Doublets
seuDoub <- lapply(seuDoub, function(seuDoub) {
    DefaultAssay(seuDoub) <- "SCT"
    seu = seuDoub %>%
        ## Run principal component analysis
        RunPCA(npcs = 50, verbose = F) %>%
        ## Embed in two dimensions (TSNE)
        RunTSNE(dims.use = 1:50) %>%
        ## Embed in two dimensions (UMAP)
        RunUMAP(dims = 1:50)
    ## pK Identification 
    sweep.res.list_seu <- paramSweep_v3(seu,PCs = 1:50,sct=TRUE)
    sweep.stats_seu <- summarizeSweep(sweep.res.list_seu, GT = FALSE)
    bcmvn <- find.pK(sweep.stats_seu)
    pK  <- bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]
    pK <- as.numeric(as.character(pK))
    nExp_poi <- round(0.039*nrow(seu@meta.data))  
    ## Run DoubletFinder with varying classification stringencies
    doubletFinder_v3(seu, PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi,
                     reuse.pANN = FALSE,sct=TRUE)
})

## Add Doublet information to Seurat Object
doubs <- lapply(seuDoub, function(x) {
    x  <- x@meta.data[10]
    colnames(x)  <-  "Doublets"
    x$cell  <-  rownames(x)
    x
})

doubs <- do.call(rbind,doubs)

seu@meta.data$Doublets <- doubs$Doublets[match(doubs$cell, rownames(seu@meta.data))]

seu <- subset(seu,
              Doublets == "Singlet")             

################################################################################
##  Filter, Normalize and process Seurat Object
###############################################################################
## Filter cells by numbers of genes and % mitochondria
seu <- subset(x = seu,
              subset = nFeature_RNA > min.genes & nFeature_RNA < max.genes &
                  nCount_RNA < max.UMI & nCount_RNA > min.UMI &
                  percent.mt < max.mito)

p2 <- VlnPlot(object = seu,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by= "orig.ident",
              pt.size = 0.25,
              cols=cols,
              ncol=3)

## Quantify genes and UMI counts
genePlot <- ggplot(seu@meta.data, aes(x=nFeature_RNA)) +
    geom_histogram(binwidth=10) +
    labs(title = "#Cells by expressed Genes - raw values")

umiPlot <- ggplot(seu@meta.data, aes(x=nCount_RNA)) +
    geom_histogram(binwidth=10) +
    labs(title = "#UMI counts - raw values")

seu <- SCTransform(object = seu,
                   variable.features.n = 3000)

## Quantify genes and UMI counts
geneNormPlot <- ggplot(seu@meta.data, aes(x=nFeature_SCT)) +
    geom_histogram(binwidth=10) +
    labs(title = "#Cells by expressed Genes - SCTransform")

umiNormPlot <- ggplot(seu@meta.data, aes(x=nCount_SCT)) +
    geom_histogram(binwidth=10) +
    labs(title = "#UMI counts - SCTransform")

## Save Plots
pdf("./2_SeuratPub/figures/SeuObj_FilterNomrQuant.pdf",width=15,height=7)
plot(p1)
plot(p2)
plot_grid(genePlot, umiPlot, align= "h")
plot_grid(geneNormPlot, umiNormPlot, align= "h")
dev.off()

## determine # of dimensions
set.seed(775)
seu <- RunPCA(object = seu, verbose = FALSE) %>%
    RunTSNE(object = seu, verbose = FALSE) %>%
    RunUMAP(object = seu, dims = 1:50,verbose=FALSE)%>%
    FindNeighbors(object = seu, dims = 1:50,verbose=FALSE) %>%
    FindClusters(seu, verbose=FALSE)

################################################################################
## Define cell type markers and print information for non-R analysis
################################################################################
## Remove small clusters
clustSize <- table(seu@meta.data$seurat_clusters)
nonSmallClusters <- names(clustSize[clustSize > 50]) ##Remove clusts w/ >50 nuc
seu <- subset(seu, idents = nonSmallClusters)
seu <- droplevels(seu)

## Get SCT expression of marker genes/cluster
DefaultAssay(seu) <- "SCT"
markers <- c("Rbfox3", ## Neuron
             "Slc17a7", ## Excitatory
             "Gad1", "Npy","Dlx1as", ## Inhibitory
             "Vip", "Sst", "Pvalb", ## Inhibitory
             "Olig1", "Ppp1r14a","Cldn11","Tspan2", ## Oligodendrocytes
             "Myt1","Pdgfra","Cspg4", ## Oligo-precursors
             "Aldoc", "Slc7a10", "Grin2c", "Fgfr3", ## Astrocytes
             "Cx3cr1", "Fcrls", "Ptprc", "Trem2", ## Microglia
             "Prox1","Gpr161","Il16") ## Locations within the hippocampus


dpdata <- DotPlot(seu, features=markers, col.min=0)
dpdata <- dpdata$data
dpdata_id <- split(dpdata,dpdata$id)
mybool <- do.call(cbind,lapply(dpdata_id,function(x) { x$avg.exp.scaled }))
clusters <- data.frame(mybool)
clusters$genes <- dpdata_id[[1]]$features.plot

## output data frame of marker expression
## writexl::write_xlsx(clusters,
##                     path="./2_SeuratPub/data/Cluster_ExpressionSCT.xlsx")

################################################################################
## Assign Clusters
## At this point, I manually assigned cell-types based on their combinations of
## markers in the excel file.  The marker expression data can be found in 
## ./files/Cluster_Expression.xlsx in the "Marker Expression" sheet and the 
## cluster assignments based on that expression is in the "ClusterInfo" sheet.
################################################################################
## Load Cluster Identities Assigned Manually
clusterNames <- readxl::read_excel("./files/Cluster_Expression.xlsx",
                                   sheet = "ClusterInfo")

## Reassign cluster ids to pruned ids
Idents(seu) <- seu@meta.data$seurat_clusters
pruned.ids <- clusterNames$Pruned_labels
names(x = pruned.ids) <- levels(x = seu)
seu <- RenameIdents(object = seu, pruned.ids)
seu[["cluster"]] = Idents(seu)

## Order clusteres
Idents(seu) <-
    factor(x = Idents(seu), levels = c("Excitatory Neuron - DG",
                                       "Excitatory Neuron - CA1",
                                       "Excitatory Neuron - CA3",
                                       "Excitatory Neuron",
                                       "Inhibitory Neuron",
                                       "Oligodendrocyte",
                                       "Oligo-precursor",
                                       "Astrocyte",
                                       "Microglia",
                                       "NA"))

## Save results
saveRDS(seu, "./files/SeuratObject.rds")
