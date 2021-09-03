## Load required libraries
library(Seurat)

## Load Data
seu <- readRDS("./2_SeuratPub/data/SeuratObject.rds")
DefaultAssay(seu) <- "RNA"

## Load Colors
cols <- pal_jco()(10)
names(cols) <- levels(unique(Idents(seu)))

################################################################################
## Run DE for each cell type
################################################################################
## Run wilcox DE analysis for each cluster in each seurat object
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
saveRDS(lrt, "./2_seuratPub/data/DE_logfc0_LR.rds")
