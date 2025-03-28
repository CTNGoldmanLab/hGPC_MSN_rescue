########################################################################
########################################################################
### convert_mouse_to_human
########################################################################
######################################################################### 
# Function to return homologous genes 
# Variables: 
#--gene_list              vector, list of genes 
#--mouse_human_genes      mouse_human_genes database
convert_mouse_to_human <- function(gene_list, mouse_human_genes){
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  
  return (output)
} # End convert_mouse_to_human

########################################################################
########################################################################
### convert_HUMAN_to_mouse
########################################################################
######################################################################### 
# Function to return homologous genes 
# Variables: 
#--gene_list              vector, list of genes 
#--mouse_human_genes      mouse_human_genes database - hard-coded and will have to download this
convert_HUMAN_to_mouse <- function(gene_list, mouse_human_genes){
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="human"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      mouse_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(mouse_gene in mouse_genes){
        output = append(output,mouse_gene)
      }
    }
  }
  
  return (output)
} # End convert_HUMAN_to_mouse
########################################################################
########################################################################
### relax_geneSet.f
########################################################################
######################################################################### 
# Function to return a data frame of each gene and the term it's mapped to 
# Variables: 
#--gsea_df              df, results from fgsea or egos
#--geneset_colnames     string, column where genes are arranged into a vector or a list 
#--term_colnames        string, column where the term is kept 
#--convert              string, hard-coded function, choice of None, toHuman, toMouse, if any of the convert function should be use 
#--delimiter            choice of ";" for fgsea result or "\\/" for egos results
relax_geneSet.f <- function(gsea_df, geneset_colnames, term_colnames, convert = "None", delimiter = ";") {
  syngo_pIDs <- gsea_df[,term_colnames, drop = TRUE]
  pway_relaxed <- data.frame(Gene = character(), Term = character())
  
  for (pway in syngo_pIDs) {
    genes <- gsea_df[which(gsea_df[, term_colnames,drop = TRUE] == pway), geneset_colnames, drop = TRUE]
    # split into vectors 
    genes <- strsplit(genes, delimiter)[[1]]
    if (convert == "toHuman") {
      genes <- convert_mouse_to_human(genes, mouse_human_genes)
    } else if (convert == "toMouse") {
      # Turn into mouse gene symbol (Remember that we had to switch to human for SynGO)
      genes <- convert_HUMAN_to_mouse(genes, mouse_human_genes)
    }
    
    # make data frame with geneName and pway 
    tmp <- data.frame(Gene = genes, Term = rep(pway, length(genes)))
    pway_relaxed <- rbind(pway_relaxed, tmp)
  }
  return(pway_relaxed)
} # END relax_geneSet.f
########################################################################
########################################################################
### get_diff.f
########################################################################
######################################################################### 
# Wrapper function to extract diff expression and subsequent clusterProfiler 
# Variables: 
#--object              seurat object, to test diff on - use MAST by default 
#--ident1,ident2       string, for Seurat FindMarkers 
#--fc_cut, padj_cut    numeric, cut-off to get signif dysregulated genes 
#--
get_diff.f <- function(object, ident1, ident2, fc_cut, padj_cut) {
  # Get everything 
  res <- FindMarkers(object, ident.1 = ident1, ident.2 = ident2, test.use = "MAST", verbose = FALSE, logfc.threshold = 0)
  #
  egosU <- enrichGO(gene = row.names(res)[res$p_val_adj < padj_cut & res$avg_log2FC > fc_cut], keyType = "SYMBOL",
                    OrgDb = org.Mm.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 1e-4,
                    qvalueCutoff = 0.01)
  egosD <- enrichGO(gene = row.names(res)[res$p_val_adj < padj_cut & res$avg_log2FC < -fc_cut], keyType = "SYMBOL",
                    OrgDb = org.Mm.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 1e-4,
                    qvalueCutoff = 0.01)
  return(list(MAST=res, egos_up = egosU, egos_dn = egosD))
} # END get_diff.f

########################################################################
########################################################################
### make_geneAnno.f
########################################################################
########################################################################
# Wrapper function to to get annotation.df for heatmaps 
# Variables: 
#--syngo_pIDs               vector, term names (fx: synapse organization)
#--pway_synGO.df            df, relaxed df from synGo results 
#--synGO.df                 df, synGO results
#--pathway_to_choose_from   c(H, I, G) in this specific case 
#--color_to_choose_from     vector of color 
#   > head(synGO)
# # A tibble: 6 × 13
#   `GO term ID` `user interface reference…` `GO domain` `GO term name` `GO term name …` `GO parent ter…` `GSEA p-value` `GSEA 'gene cl…`
#   <chr>        <chr>                       <chr>       <chr>          <chr>            <chr>                     <dbl>            <dbl>
# 1 GO:0099536   H1                          BP          synaptic sign… ├─ synaptic sig… SYNGO:synprocess       1.49e-14         2.24e-13
make_geneAnno.f <- function(syngo_pIDs, pway_synGO.df, synGO.df, pathway_to_choose_from, color_to_choose_from) {
  # Map to parent terms 
  x <- as.data.frame(table(pway_synGO.df$Term))
  x$reference <- plyr::mapvalues(x$Var1, synGO.df$`GO term name`, synGO.df$`user interface reference code`, warn_missing = FALSE) %>% substring(1, 1)
  pway_synGO.df2 <- pway_synGO.df
  pway_synGO.df2$reference <- plyr::mapvalues(pway_synGO.df2$Term, x$Var1, x$reference) 
  pway_synGO.df2 <- pway_synGO.df2[, c("Gene", "reference")][!duplicated(pway_synGO.df2[, c("Gene", "reference")]),]
  # Asign colors for the unique genes 
  # For non-unique genes, assign "black"
  tmp <- pway_synGO.df2
  tmp <- table(tmp) 
  tmp <- apply(tmp, 1, sum) %>% as.data.frame()
  
  annot_row.df <- data.frame(row.names = character(), Term = character(), Color = character())
  for (g in row.names(tmp)) {
    if (tmp[g,1] > 1) {
      token <- "multiple"
      token_color <- "black"
    } else {
      token_df <- pway_synGO.df2
      token <- token_df$reference[which(token_df$Gene == g)]
      token_color <- plyr::mapvalues(token, pathway_to_choose_from, color_to_choose_from, warn_missing = FALSE)
    }
    annot_row.df <- rbind(annot_row.df, data.frame(row.names = g, Term = token, Color = token_color))
  }
  return(annot_row.df)
} # END make_geneAnno.f
########################################################################
########################################################################
### lmp
########################################################################
########################################################################
# Function to extract p value from linear model object
# Source: http://www.gettinggeneticsdone.com/2011/01/rstats-function-for-extracting-f-test-p.html
# STAR
lmp <- function(modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
} # END lmp

########################################################################
########################################################################
### calculate_relative_importance
########################################################################
# Function to calculate relative importance of known covariates in the module
# eigengene (ME) variation. The importance is measured by the relative contribution
# of each covariate towards the explanatory value of the additive R^2. The results
# are saved in output files. See documentation of relaimpo package.
# Function returns a data frame of each covariate by ME - and their confidence intervals 
#
# Parameters:
#  suffix (string):     suffix saved results 
#  fit (lm):            linear model object containing the additive model with ME as
#                       the response variable and known covariates as predictors
#  cov_keep (vector):   a vector of covariate considered for this linear regression. hard-coded.
# STAR
calculate_relative_importance <- function(fit, suffix) {
  
  bootfit <- NULL
  # may come back as NULL when dependent covariates produce a singlularity error when
  # calculating the bootstrap intervals
  while (is.null(bootfit)) {
    bootfit <- boot.relimp(fit, b = 200, type = c("lmg", "first", "last"), rank = TRUE, diff = TRUE, rela = TRUE, fixed = TRUE)
  } # end while
  
  # plot the bootfit results
  # pdf(paste0("./plots/Relative Covariate Importance/", suffix, "_Booteval Plot.pdf"),
  #     height = 6, width = 8, useDingbats = FALSE)
  # plot(booteval.relimp(bootfit))
  # graphics.off()
  
  bootout <- booteval.relimp(bootfit)
  
  res <- rbind(t(rbind(bootout@lmg, bootout@lmg.lower, bootout@lmg.upper)),
               t(rbind(bootout@first, bootout@first.lower, bootout@first.upper)), 
               t(rbind(bootout@last, bootout@last.lower, bootout@last.upper)))
  res <- as.data.frame(res)
  colnames(res)  <- c("confidence", "lower", "upper") 
  res$covariate  <- paste(rep(cov_keep, 3), rep(c("lmg", "first", "last"), each = length(cov_keep)), sep = "_")
  res$TF         <- rep(suffix, nrow(res))
  
  # save output to text files
  # write.table(res,
  #             paste0("./tables/Relative Covariate Importance/", suffix, ".txt"),
  #             sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # extract data frame to plot covariate contribution to ME
  res1           <- res[1:length(cov_keep),] # only interested in the lmg method
  res1$covariate <- gsub("_lmg", "", res1$covariate)
  return(res1)
} # end calculate_relative_importance()