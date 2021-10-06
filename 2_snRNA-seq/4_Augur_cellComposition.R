## #############################################################################
## Date:        April 2021
## Author:      Allison M. Burns
## Filename:    /2_snRNA-seq/4_Augur_cellComposition.R
## Project:     Epigenetic Priming - snRNA-seq analysis
## Description: Code for running Augur and the permutation test.
##              Augur (prioritization of a population’s responsiveness to
##              experimental perturbation) tests to determine the extent to
##              which HDACi induces differential expression in different cell
##              types, without being biased by numbers of nuclei.
##              And the permutation test determines HDACi alters the numbers
##              of nuclei within a give cell-type.
## #############################################################################

library(Seurat)
library(Augur)
library(scProportionTest)
library(ggplot2)

################################################################################
## Run Augur
################################################################################
## Load datasets
seu <- readRDS("./files/SeuratObject.rds")
DefaultAssay(seu)  <-  "RNA"
unique(Idents(seu))

## Run Augur on assigned identities
aug <-  calculate_auc(seu, label_col = "drug",cell_type_col = "cluster")
saveRDS(aug, "./files/Augur_seurat_clusters.rds")

## Internal augur lollipops
plot_lollipop(aug) +
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10))

################################################################################
## Run scProportions to determine changes in cell composition
################################################################################
## scProportions test - https://www.biostars.org/p/457243/
prop_test <- sc_utils(seu)

prop_test <- permutation_test(
    prop_test, cluster_identity = "cluster",
    sample_1 = "Veh", sample_2 = "HDACi",
    sample_identity = "drug"
)


permutation_plot(prop_test, log2FD_threshold = log2(1.5)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(family = "Helvetica",color = "black", size = 12),
          axis.title = element_text(family = "Helvetica",color = "black", size = 12),
          axis.title.y = element_blank()
          )

