# Use this script to export matrix for pySCENIC
integrated_object <- readRDS("outs/RDS/integrated_object.RDS")
##### Export matrix for SCENIC
scenic <- integrated_object[["RNA"]]@counts
dim(scenic). # 34039 53259
# Filter to retain only genes that are expressed more than minimal counts per gene threshold
minCountsPerGene <- 3*0.01*ncol(scenic)
sum <- apply(scenic, 1, sum)
keep <- names(sum)[which(sum > minCountsPerGene)]
# Filter to retain only genes that are expressed in more than the minimal number of cells
minSamples <- ncol(scenic)*0.01
zeros <- apply(scenic, 1, function(x) sum(x>0))
keep2 <- names(zeros)[which(zeros > minSamples)]
keep <- keep[keep %in% keep2]
scenic <- as.matrix(scenic[keep,])
dim(scenic) # 13682 53259

write.csv(scenic, paste0("outs/tables/forScenic.csv"), quote = FALSE)
