## Load required libraries
library(Seurat)
library(ggsci)
library(ggplot2)
library(ggrastr)

## Load Data
seu <- readRDS("./2_SeuratPub/data/SeuratObject.rds")
DefaultAssay(seu) <- "RNA"
lrt <- readRDS("./2_seuratPub/data/DE_logfc0_LR.rds")

## Load Colors
cols <- pal_jco()(10)
names(cols) <- levels(unique(Idents(seu)))

## Define Differential Enrichment
fc <- 1
pval <- 0.05

################################################################################
##  Remove gene groups
################################################################################
set.seed(775)
allGenes <- rownames(seu@assays$RNA@data)
up <- lapply(lrt, function(x) {x[x$p_val_adj <= pval & x$avg_log2FC >= fc,]})
dw <- lapply(lrt, function(x) {x[x$p_val_adj <= pval & x$avg_log2FC <= -fc,]})
all <- lapply(lrt, function(x) {x[x$p_val_adj <= pval & abs(x$avg_log2FC) >= fc,]})
ran <- lapply(lrt, function(x) { x[sample(seq(nrow(x)),nrow(x)/2),] })

## UMAP loop for all given gene lists
runUMAPremoval <- function(i,geneLists,seuData,fl) {
    ## Set up list information
    myList <- geneLists[[i]]
    listName <- names(geneLists[i])
    cat(paste("Filtering for seurat object for", listName ,"\n",sep=" "))
    
    ## Setup Seurat Object
    set.seed(775)
    seuData <- subset(x = seu,
                      features = allGenes[!allGenes %in% rownames(myList)])
    DefaultAssay(seuData)  <- 'SCT'
    
    ## Dimension Reduction
    seuData <- RunPCA(object = seuData, verbose = FALSE)
    seuData <- RunUMAP(object = seuData, dims = 1:50,verbose=FALSE)
    
    ## Make Figures
    cols <- c("#8B829C", "#5A486F","#99CCBA", "#689D6C") ## telegraph
    
    ## Make Figures
    dp1 <- rasterize(DimPlot(object = seuData,
                             reduction = "umap",
                             group.by = "drug",
                             cols=cols[c(2,4)],
                             label = FALSE),
                     dpi = 400)
    
    dp2 <- DimPlot(object = seuData,
                   reduction = "umap",
                   label = TRUE,
                   label.size = 3.5,
                   pt.size = 0.5,
                   order = TRUE) +
        NoLegend()
    dp2$layers[[2]]$aes_params$fontface <- "bold"
    
    dp1$layers[[2]] <- dp2$layers[[2]]
    
    pdf(paste("./2_SeuratPub/figures/DivideConquer/",fl,"/",listName,".pdf",sep = ""),
        width = 5, height = 5)
    print(
        dp1 +
        labs(title = paste("UMAP after removing up-regulated genes\nfrom",
                           listName ,"and re-clustering",sep = " ")) +
        NoLegend() +
        theme (plot.title = element_text(hjust = 0.5,size = 12, face = "bold"),
               plot.caption = element_text(face = "italic"),
               ##plot.title = element_blank(),
               text = element_text(family = "Helvetica",face = "bold"),
               legend.position = c(0.05, 0.9))
    )
    dev.off()
}

upUMAP <- lapply(seq(length(up)), runUMAPremoval, up, seu, "up_genes")
dwUMAP <- lapply(seq(length(dw)), runUMAPremoval, dw, seu, "dw_genes")
ranUMAP <- lapply(seq(length(ran)), runUMAPremoval, ran, seu, "ran_genes")
allUMAP <- lapply(seq(length(all)), runUMAPremoval, all, seu, "all_DEGs")
