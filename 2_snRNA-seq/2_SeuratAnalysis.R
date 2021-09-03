## Load required libraries
library(Seurat)
library(dplyr)
library(sctransform)
library(DoubletFinder)
library(ggplot2)
library(cowplot)

################################################################################
## Load count matrices for all files
################################################################################
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

## Set Plotting Colors
cols <- c("#8B829C", "#5A486F","#99CCBA", "#689D6C") ## telegraph

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
seu <- RunPCA(object = seu, verbose = FALSE) %>%
    RunTSNE(object = seu, verbose = FALSE) %>%
    RunUMAP(object = seu, dims = 1:50,verbose=FALSE)%>%
    FindNeighbors(object = seu, dims = 1:50,verbose=FALSE) %>%
    FindClusters(seu, verbose=FALSE)

################################################################################
## Name clusters
################################################################################
## Load Cluster Identities Assigned Manually
clusterNames <- readxl::read_excel("./2_SeuratPub/data/Cluster_Expression.xlsx",
                                   sheet = "ClusterInfo")

## Save Full cell type information
full.ids <- paste(clusterNames$Cell_type,"-",clusterNames$Cell_location,sep=" ")
names(x = full.ids) <- levels(x = seu)
seu <- RenameIdents(object = seu, full.ids)
seu[["full_clusterID"]] = Idents(seu)

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
saveRDS(seu, "./2_SeuratPub/data/SeuratObject.rds")
