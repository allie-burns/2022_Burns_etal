## #############################################################################
## Date:        April 2021
## Author:      Allison M. Burns
## Filename:    /2_snRNA-seq/3_DEanalysis.R
## Project:     Epigenetic Priming - snRNA-seq analysis
## Description: Differential expression analysis between HDACi and Veh treated
##              samples for each cell-type in the Seurat object. This outputs the
##              log2FC and p-value for all expressed genes (in > 10% of cells)
##              within each cell type. 
## #############################################################################

library(Seurat)

## Load Data
seu <- readRDS("./files/SeuratObject.rds")
DefaultAssay(seu) <- "RNA"

################################################################################
## Run DE for each cell type
################################################################################
## Run logistic regression framework DE analysis for each cluster 
clusters <- levels(Idents(seu))

DEbyCluster <- function(clusters,seuData) {
    seu.ident <- subset(x = seuData, idents = clusters)
    Idents(object = seu.ident) <- "drug"
    
    FindMarkers(seu.ident,
                ident.1 = 'HDACi',
                ident.2 = 'Veh',
                logfc.threshold=0,
                test.use = "LR", ## Comment out for wilcox
                latent.vars = "day", ## Comment out for wilcox
                min.pct=0.1) ## gene is considered expressed if in >= 10% of cells
}

lrt <- lapply(clusters,DEbyCluster,seu)
names(lrt)  <- clusters
saveRDS(lrt, "./files/DE_logfc0_LR.rds")
