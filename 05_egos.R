library(clusterProfiler)
library(org.Mm.eg.db)

x_msn <- readRDS("outs/RDS/DEG_MSN.RDS")
# Check that the down-regulated genes are not in the up-regulated list as well 
y1 <- unlist(x_msn[c("d1_dn", "d2_dn")]) %>% unique()
y2 <- unlist(x_msn[c("d1_up", "d2_up")]) %>% unique()
table(y2 %in% y1); table(y1 %in% y2)


egos_MSN_U <- enrichGO(gene = y2, keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 1e-4,
                       qvalueCutoff = 0.01)
egos_MSN_D <- enrichGO(gene = y1, keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 1e-4,
                       qvalueCutoff = 0.01)
save(egos_MSN_U, egos_MSN_D, file = "outs/RDS/MSN_union_egos.RData")