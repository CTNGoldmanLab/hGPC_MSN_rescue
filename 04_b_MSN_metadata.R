library(Seurat)
library(SeuratDisk)

Convert("outs/h5ad/msn.h5ad", "h5seurat", assay="RNA")
integrated_objectMSN <- LoadH5Seurat("outs/h5ad/msn.h5seurat", assays = "RNA")
# SCENIC results 
# FOR THE MSN ALONE
auc <- read.csv("pySCENIC_output_from_03_c/auc.txt", row.names=1)
colnames(auc) <- gsub("\\..*", "", colnames(auc))
aucMSN <- auc[row.names(auc) %in% colnames(integrated_objectMSN),]
# Have the scale.data with AUC (AUC is scaled from 0 to 1) & find clusters to overlay with UMAP
integrated_objectMSN[["AUC"]] <- CreateAssayObject(data = t(aucMSN))
integrated_objectMSN[["AUC"]]@scale.data <- integrated_objectMSN[["AUC"]]@data

#------- Assignment of healthy vs sick D1 and D2 
p1 <- DimPlot(integrated_objectMSN, group.by = "cellTypeInt")
p2 <- DimPlot(integrated_objectMSN, group.by = "mouse_engraft", cols = c("blue", "navy", "grey80", "grey30"), split.by = "mouse")
p3 <- DimPlot(integrated_objectMSN, group.by = "leiden", cols = c("firebrick4", "goldenrod4"))
cowplot::plot_grid(p2, cowplot::plot_grid(p1, p3, ncol = 2), ncol = 1)
# 0 is the "healthy cluster", and 1 is the "sick" cluster 
msnStatus <- c()
for (i in row.names(integrated_objectMSN@meta.data)) {
  subCell_token <- integrated_objectMSN$cellTypeInt[i]
  subLeiden_token <- integrated_objectMSN$leiden[i]
  if (subCell_token == "D1_MSN") {
    if (subLeiden_token == "0") {
      msnStatus[i] <- "D1_WT"
    } else if (subLeiden_token == "1") {
      msnStatus[i] <- "D1_HD"
    }
  } else if (subCell_token == "D2_MSN") {
    if (subLeiden_token == "0") {
      msnStatus[i] <- "D2_WT"
    } else if (subLeiden_token == "1") {
      msnStatus[i] <- "D2_HD"
    }
  }
}
all.equal(row.names(integrated_objectMSN@meta.data), names(msnStatus))
integrated_objectMSN$msnStatus <- msnStatus
#------- RESTORE VS SICK 
# By D1 and D2, find out what are the genes that are differentially expressed between R62+GPC that are spending their time in the WT clusters, compared to the rest 
met <- integrated_objectMSN@meta.data
restore <- c()
for (i in 1:nrow(integrated_objectMSN@meta.data)) {
  xme <- met$mouse_engraft[i]
  xct <- met$msnStatus[i] 
  if (xct == "D1_WT" & xme == "R62_hGPC") { # if engrafted and spending time with WT cells 
    token = "D1_restore"
  } else if (xct == "D2_WT" & xme == "R62_hGPC") {
    token = "D2_restore"
  } else if (xct == "D1_HD" & xme == "R62_hGPC") { # if engrafted and still sick 
    token = "D1_sick"
  } else if (xct == "D2_HD" & xme == "R62_hGPC") {
    token = "D2_sick"
  } else if (xct == "D1_HD" & xme == "R62_Unengrafted") { # If unengrafted and sick 
    token = "D1_sick"
  } else if (xct == "D2_HD" & xme == "R62_Unengrafted") {  
    token = "D2_sick"
  } else { # All other D*_ that are healthy - These come from all four groups, with just miniscule cells from R62_Unengrafted here 
    token = ""
  }
  restore <- c(restore, token)
}

integrated_objectMSN$restore <- restore
###---- RESTORE2 
# While R62_hGPC MSNs were split between WT-like and R62_Unengrafted-like cells, the ones that were designated as "sick" actually expressed genes in-between. So we would like to exclude these cells for downstream analysis for a cleaner list of DE genes 
# Therefore: The following identities will be defined as : 
# Baseline is R62_Unengrafted and annotated as sick cells. 
# Restore is R62_hGPC and annotated as restore. 
# We exclude the R62_hGPC that were annotated as "sick", since it seems like these guys were "in the process of getting better", and thus muddied our analysis. 
tmp1 <- ifelse(integrated_objectMSN$mouse_engraft == "R62_hGPC" & grepl("restore", integrated_objectMSN$restore), "restore", "")
tmp2 <- ifelse(integrated_objectMSN$mouse_engraft == "R62_Unengrafted" & grepl("sick", integrated_objectMSN$restore), "baseline", "")
integrated_objectMSN$restore2 <- paste0(tmp1, tmp2)
# saveRDS(integrated_objectMSN, "outs/RDS/integrated_objectMSN.RDS")
