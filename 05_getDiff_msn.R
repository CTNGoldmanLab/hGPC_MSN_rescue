# Set libPaths to global 
library(Seurat)
library(sctransform)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggthemes)
library(knitr)
library(clusterProfiler)
library(GO.db)
library(org.Mm.eg.db)
library(fgsea)
library(GeneOverlap)
library(VennDiagram)
library(jsonlite)
library(GeneOverlap)
library(glmnet)
library(relaimpo)

source("helper_data_wrangling_functions.R")
source("helper_plot_functions.R")

integrated_object <- readRDS("outs/RDS/integrated_object.RDS")
DefaultAssay(integrated_object) <- "RNA" # "RNA"

####  Neuron populations between:
##### Neuron R62 unengrafted and WT unengrafted (rw_un)
markers_msn_pair <- list()
pair_toCompare <- list(rw_un = c("R62_Unengrafted", "WT_Unengrafted"))
for (ct in c("D1_MSN", "D2_MSN")) {
  subOb <- subset(integrated_object, idents = ct)
  subOb <- SetIdent(subOb, value="mouse_engraft")
  for (ident in names(pair_toCompare)) {
    markers_msn_pair[[ident]][[ct]] <- get_diff.f(subOb, 
                                                  ident1 = pair_toCompare[[ident]][1], 
                                                  ident2 = pair_toCompare[[ident]][[2]], 
                                                  fc_cut = 0.5, padj_cut = 1e-3)
  }
}
# saveRDS(markers_msn_pair, "outs/RDS/MAST_MSN_pairwise.RDS")