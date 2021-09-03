library(DESeq2)

## Download DE results
BR  <- "hpp"
results <- readRDS(paste("./data/DE_conditioning/DEseq2",BR,"VCcontrol.rds",sep = "_"))
results <- results[1:3] ## remove HDACi comparison
names(results) <- c("VS","HC","HS")

################################################################################
## Set up FC information
################################################################################
## Get combined of all genes across all comparisons
gns <- unlist(lapply(results,rownames))
gns <- unique(gns)

## order results to match in all groups
results <- lapply(results, function(x) {
    x <- x[match(gns,rownames(x)),]
    rownames(x)  <-  gns
    x
})

##build dataframe
##set all insig values to 0 (appear as grey in heatmap)
log2FC <- lapply(results, function(x) {
    FC <- x$log2FoldChange
    names(FC)  <-  rownames(x)
    FC[is.na(FC)]  <-  0
    FC[x$padj > 0.05]  <- 0
    FC
})

log2FC$VC <- rep(0, length(log2FC[[1]]))
log2FC <- log2FC[c(4,1,2,3)]

FC <- do.call(cbind,log2FC)

################################################################################
## Define clusters
################################################################################
clustFC <- data.frame(FC)
clustFC$gene <- rownames(clustFC)

d = 0.2

clusters <- lapply(seq(1:nrow(clustFC)), function(i) {
    VC  <-  clustFC[i,1]
    VS  <-  clustFC[i,2]
    HC  <-  clustFC[i,3]
    HS  <-  clustFC[i,4]
    
    if(VS >= VC + d) {
        if(HC >= VS + d){
            if(HS >= HC + d){
                c = "1.1.1"
            }else if(HS < HC + d & HS > HC - d){
                c = "1.1.2"
            }else if(HS <= HC - d){
                c = "1.1.3"
            }
        }else if(HC < VS + d & HC > VS -d){
            if(HS >= HC + d){
                c = "1.2.1"
            }else if(HS < HC + d & HS > HC - d){
                c = "1.2.2"
            }else if(HS <= HC - d){
                c = "1.2.3"
            }
        }else if(HC <= VS - d){
            if(HS >= HC + d){
                c = "1.3.1"
            }else if(HS < HC + d & HS > HC - d){
                c = "1.3.2"
            }else if(HS <= HC - d){
                c = "1.3.3"
            }
        }
    }else if(VS < VC + d & VS > VC - d) {
        if(HC >= VS + d){
            if(HS >= HC + d){
                c = "2.1.1"
            }else if(HS < HC + d & HS > HC -d){
                c = "2.1.2"
            }else if(HS <= HC - d){
                c = "2.1.3"
            }
        }else if(HC < VS + d & HC > VS -d){
            if(HS >= HC + d){
                c = "2.2.1"
            }else if(HS < HC + d & HS > HC -d){
                c = "2.2.2"
            }else if(HS <= HC - d){
                c = "2.2.3"
            }
        }else if(HC <= VS - d){
            if(HS >= HC + d){
                c = "2.3.1"
            }else if(HS < HC + d & HS > HC -d){
                c = "2.3.2"
            }else if(HS <= HC - d){
                c = "2.3.3"
            }
        }
    }else if(VS <= VC - d) {
        if(HC >= VS + d){
            if(HS >= HC + d){
                c = "3.1.1"
            }else if(HS < HC + d & HS > HC -d){
                c = "3.1.2"
            }else if(HS <= HC - d){
                c = "3.1.3"
            }
        }else if(HC < VS + d & HC > VS -d){
            if(HS >= HC + d){
                c = "3.2.1"
            }else if(HS < HC + d & HS > HC -d){
                c = "3.2.2"
            }else if(HS <= HC - d){
                c = "3.2.3"
            }
        }else if(HC <= VS - d){
            if(HS >= HC + d){
                c = "3.3.1"
            }else if(HS < HC + d & HS > HC -d){
                c = "3.3.2"
            }else if(HS <= HC - d){
                c = "3.3.3"
            }
        }
    }
})

clustFC$cluster <- unlist(clusters)


saveRDS(object = clustFC,
        file = paste("./data/DE_conditioning/clusterFC_",BR,".rds",sep = ""))

