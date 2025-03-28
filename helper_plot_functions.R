########################################################################
########################################################################
### draw_heatmap.f
########################################################################
########################################################################
# Function to return pheatmap-like but with ggplot2
# Variables: 
#--smple_df              df, rownames of sample, and colnames of sample meta data 
#--expr_mat              matrix, rownames of gene, colnames of samples 
main_theme <- theme(axis.text = element_text(colour = "black"), axis.title = element_text(colour = "black"), panel.grid = element_blank()) 
draw_heatmap.f <- function(expr_mat, smple_df, scale = "none", showSample = FALSE, showGene = TRUE) {
  ##########---------------
  # Scale by row if option is selected 
  if (scale == "row") {
    expr_mat <- t(apply(expr_mat, 1, function(x) (x-mean(x))/sd(x)))
  }
  ##########---------------
  # Hierachical cluster 
  hcRow <- hclust(dist(expr_mat, method = "euclidean"))
  ##########---------------
  # Plot main tile 
  expr_mat.long <- expr_mat %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") 
  expr_mat.long <- tidyr::gather(expr_mat.long, "sample", "expression", 2:ncol(expr_mat.long))
  # Prevent ggplot auto-ordering 
  expr_mat.long$gene   <- factor(expr_mat.long$gene,   levels = row.names(expr_mat)[hcRow$order])
  expr_mat.long$sample <- factor(expr_mat.long$sample, levels = colnames(expr_mat))
  main.p <- ggplot(expr_mat.long, aes(x=sample, y = gene, fill = expression)) + geom_tile(color="black") + theme_minimal_grid() + main_theme + labs(y = "", x = "") + scale_y_discrete(position = "right") + theme(axis.text.y = element_text(size = 9), legend.title = element_text(size=7.5), legend.text = element_text(size = 7), legend.position = "right", axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1))
  
  if (!showSample) {
    main.p <- main.p + theme(axis.text.x = element_blank())
  }
  if (!showGene) {
    main.p <- main.p + theme(axis.text.y = element_blank())
  }
  if (scale == "row") {
    main.p <- main.p + labs(fill = "Row Z-score")
    main.p <- main.p + scale_fill_gradient2(low = "navyblue", high = "firebrick4", mid = "white", midpoint = 0)
  } else {
    main.p <- main.p + scale_fill_viridis(option = "B", begin = 0) 
  }
  return(main.p)
} # end draw_heatmap.f


########################################################################
########################################################################
### plot_egos.f
########################################################################
########################################################################
# Wrapper function to extract signif egos and plot from markers_glial 
plot_egos.f <- function(markers_list, wrap_width = 25) {
  # markers_list: list of cell type with enwrapped list of (MAST = res, egos_up = egosU, egos_dn = egosD)
  # function return a plot list of top 5 upregulated and top 5 downregulated pways 
  markers_temp <- list()
  markers_temp[["dn"]] <- lapply(markers_list, function(x) x[["egos_dn"]]@result %>% mutate(logP = -log10(p.adjust)) %>% top_n(-4, p.adjust))
  markers_temp[["up"]] <- lapply(markers_list, function(x) x[["egos_up"]]@result %>% mutate(logP = -log10(p.adjust)) %>% top_n(-4, p.adjust))
  
  # For each cell type, go through and extract data to plot (dtp)
  plist <- list()
  for (i in names(markers_temp[["up"]])) {
    dtp <- rbind(markers_temp[["dn"]][[i]], markers_temp[["up"]][[i]])
    # Factor for terms 
    rearrange_x <- c(markers_temp[["dn"]][[i]]$Description[order(markers_temp[["dn"]][[i]]$p.adjust, decreasing = TRUE)],
                     markers_temp[["up"]][[i]]$Description[order(markers_temp[["up"]][[i]]$p.adjust, decreasing = TRUE)])
    dtp$Description <- factor(dtp$Description, levels = rearrange_x)
    
    p <- ggplot(dtp, aes(x=Description, y = logP)) + geom_col(fill = ifelse(dtp$Description %in% markers_temp[["dn"]][[i]]$Description, "navy", "firebrick4")) + coord_flip() + theme_linedraw() 
    p <- p + labs(x = "", y = "-log10(padj)", subtitle = i) + scale_x_discrete(labels = function(x) str_wrap(x, width = wrap_width)) + theme(axis.text.y = element_text(size = 10), panel.grid = element_blank())
    plist[[i]] <- p
  }
  return(plist)
}

########################################################################
########################################################################
### plot_egos2.f
########################################################################
########################################################################
# Wrapper function to extract signif egos and plot from markers_glial 
# Have an if statement to skip egos if there's no signif
plot_egos2.f <- function(markers_list) {
  # markers_list: list of cell type with enwrapped list of (MAST = res, egos_up = egosU, egos_dn = egosD)
  # function return a plot list of top 5 upregulated and top 5 downregulated pways 
  # Inspec markers_list and make lean (i.e. remove NULL egos results since there's no signif enrichment)
  markers_list_up <- lapply(markers_list, function(x) x[["egos_up"]])
  markers_list_dn <- lapply(markers_list, function(x) x[["egos_dn"]])
  # The cell type to exclude because NULL
  exclude_up <- names(markers_list)[unlist(lapply(markers_list_up, is.null))]
  markers_list_up <- markers_list_up[!names(markers_list_up) %in% exclude_up]
  #
  exclude_dn <- names(markers_list)[unlist(lapply(markers_list_dn, is.null))]
  markers_list_dn <- markers_list_dn[!names(markers_list_dn) %in% exclude_dn]
  #
  markers_temp <- list()
  markers_temp[["dn"]] <- lapply(markers_list_dn, function(x) x@result %>% mutate(logP = -log10(p.adjust)) %>% top_n(-3, p.adjust))
  markers_temp[["up"]] <- lapply(markers_list_up, function(x) x@result %>% mutate(logP = -log10(p.adjust)) %>% top_n(-3, p.adjust))
  
  # For each cell type, go through and extract data to plot (dtp)
  plist <- list()
  for (i in names(markers_list)) {
    dtp <- rbind(markers_temp[["dn"]][[i]], markers_temp[["up"]][[i]])
    # Factor for terms
    if (is.null(markers_list[[i]]$egos_up)) {
      rearrange_x <- c(markers_temp[["dn"]][[i]]$Description[order(markers_temp[["dn"]][[i]]$p.adjust, decreasing = TRUE)])
    } else if (is.null(markers_list[[i]]$egos_dn)) {
      rearrange_x <- c(markers_temp[["up"]][[i]]$Description[order(markers_temp[["up"]][[i]]$p.adjust, decreasing = TRUE)])
    } else {
      rearrange_x <- c(markers_temp[["dn"]][[i]]$Description[order(markers_temp[["dn"]][[i]]$p.adjust, decreasing = TRUE)],
                       markers_temp[["up"]][[i]]$Description[order(markers_temp[["up"]][[i]]$p.adjust, decreasing = TRUE)])
    }
    
    dtp$Description <- factor(dtp$Description, levels = rearrange_x)
    
    p <- ggplot(dtp, aes(x=Description, y = logP)) + geom_col(fill = ifelse(dtp$Description %in% markers_temp[["dn"]][[i]]$Description, "navy", "firebrick4")) + coord_flip() + theme_linedraw() 
    p <- p + labs(x = "", y = "-log10(padj)", subtitle = i) + scale_x_discrete(labels = function(x) str_wrap(x, width = 25)) + theme(axis.text.y = element_text(size = 10), panel.grid = element_blank())
    plist[[i]] <- p
  }
  return(plist)
} # end plot_egos2.f

########################################################################
########################################################################
### draw_heatmap_v2
########################################################################
########################################################################
#####################################################################################
#####################################################################################
# draw_heatmap_v2()
# Function to return a heatmap based on the matrix provided
# Function relates to the egos gene ontology term for top_node
# Parameters:
#--annot_row             df, rownames of genes, colnames of pway/TF, with "" or exact replica of pway/TF in the values
#--smple_df              df, rownames of sample, and colnames of sample meta data 
#--expr_mat              matrix, rownames of gene, colnames of samples 
#--subtext               string, what to call the plot
draw_heatmap_v2 <- function(expr_mat, annot_row, smple_df, scale = "none", showSample = FALSE, showGene = TRUE, subtext) {
  ##########---------------
  # Scale by row if option is selected 
  if (scale == "row") {
    expr_mat <- t(apply(expr_mat, 1, function(x) (x-mean(x))/sd(x)))
  }
  ##########---------------
  # Hierachical cluster 
  hcRow <- hclust(dist(expr_mat, method = "euclidean"))
  ##########---------------
  # Plot main tile 
  expr_mat.long <- expr_mat %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") 
  expr_mat.long <- tidyr::gather(expr_mat.long, "sample", "expression", 2:ncol(expr_mat.long))
  # Prevent ggplot auto-ordering 
  expr_mat.long$gene   <- factor(expr_mat.long$gene,   levels = row.names(expr_mat)[hcRow$order])
  expr_mat.long$sample <- factor(expr_mat.long$sample, levels = colnames(expr_mat))
  main.p <- ggplot(expr_mat.long, aes(x=sample, y = gene, fill = expression)) + geom_tile(color="black") + theme_minimal_grid() + main_theme + labs(y = "", x = "") + scale_y_discrete(position = "right")
  main.p <- main.p + theme(axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1), 
                           axis.text.y = element_text(size = 7.2), 
                           legend.title = element_text(size=7.5), 
                           legend.text = element_text(size = 7), 
                           legend.position = "right", 
                           plot.margin = unit(c(2,0,0,0), "mm"))
  
  if (!showSample) {
    main.p <- main.p + theme(axis.text.x = element_blank())
  }
  if (!showGene) {
    main.p <- main.p + theme(axis.text.y = element_blank())
  }
  if (scale == "row") {
    main.p <- main.p + labs(fill = "Row Z-score", subtitle = subtext)
    main.p <- main.p + scale_fill_gradient2(low = "navyblue", high = "firebrick4", mid = "white", midpoint = 0)
  } else {
    main.p <- main.p + scale_fill_viridis(option = "B", begin = 0) + labs(fill = "Scaled\nExpression", subtitle = subtext)
  }
  ##########---------------
  # annotation
  my_color <- c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(8, "Accent"))
  # values of colors so that it stays consistent among plots
  my_color.sub <- c("ghostwhite", my_color[1:ncol(annot_row)])
  names(my_color.sub) <- c("",colnames(annot_row)[order(colnames(annot_row))])
  annot <- annot_row[which(row.names(annot_row) %in% row.names(expr_mat)),] %>% tibble::rownames_to_column(var="gene")
  annot <- tidyr::gather(annot, GO, annotation, 2:ncol(annot))
  annot$Target <- factor(annot$gene, levels = row.names(expr_mat)[hcRow$order])
  x <- unique(annot$annotation)
  x <- x[order(x)]
  annot$annotation <- factor(annot$annotation, levels = x)
  my_color.sub <- my_color.sub[x]
  h <- ggplot(annot, aes(x = GO, y = Target, fill = annotation)) + geom_tile(color = "black") + scale_fill_manual(values = my_color.sub, na.value="white") + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + theme_minimal()
  h <- h + labs(x="", y="") + theme(axis.text.y = element_blank(), 
                                    axis.text.x = element_text(color = "black", angle = 45, hjust = 1,size=8), 
                                    axis.ticks.x = element_blank(), 
                                    axis.ticks.y.left = element_blank(), 
                                    axis.ticks.y.right = element_line(color = "black"), 
                                    legend.position = "none", plot.margin = unit(c(2,0,0,0), "mm"))
  #
  return(list(h, main.p))
} # end draw_heatmap_v2
########################################################################
########################################################################
### do_heatmapSeurat.f
########################################################################
########################################################################
# Function to return DoHeatMap-like but with ggplot2
# Function plot a color bar for cell column, 
# and utilize different text color for row
# Variables: 
#--object                SeuratObject. Object will be subset based on active.ident
#--annot_df              df, gene as row names, and colnames with annotation terms and accompanying colors (hardcoded: Term, Color)
#--smple_df              df, rownames of cell, and colnames of seurat meta data 
#--features              vector, a list of genes of interest 
#--color_r               vector, carrying color indication for sample by group.by c(R62="black", WT="snow")
#--prop_size             numeric, proportion to draw samples from
do_heatmapSeurat.f <- function(object, annot_df, smple_df, features, assay = "RNA", scale = TRUE, showSample = FALSE, showGene = TRUE, color_list, prop_size=0.1) {
  ##########---------------
  # Default is to use the scaled data plot. If not, use normalized slot
  # Default is to use RNA assay. Can substitute with other options. 
  if (scale) {
    expr_mat <- object[[assay]]@scale.data 
  } else {
    expr_mat <- object[[assay]]@data
  }
  # Select by features 
  expr_mat <- expr_mat[row.names(expr_mat) %in% features,]
  # Randomly take 10% per ident 
  # What consitute an ident is dependent on smple_df
  set.seed(2208)
  tmp_cols <- colnames(smple_df)
  if (length(tmp_cols) > 1) {
    smple_df$tmp <- apply(smple_df[, tmp_cols], 1, paste, collapse = "_")
  } else {
    smple_df$tmp <- smple_df[,1]
  }
  cells_to_subset <- c()
  for (group in unique(smple_df$tmp)[order(unique(smple_df$tmp))]) {
    tmp <- sample(x=row.names(smple_df)[which(smple_df$tmp == group)], size = prop_size*sum(smple_df$tmp == group))
    cells_to_subset <- c(cells_to_subset, tmp)
  }
  # Filter to remain subset cells 
  expr_mat <- expr_mat[,cells_to_subset]
  ##########---------------
  # Hierachical cluster 
  hcRow <- hclust(dist(expr_mat, method = "euclidean"))
  ##########---------------
  # Plot main tile 
  expr_mat.long <- expr_mat %>% as.data.frame %>% tibble::rownames_to_column(var = "gene") 
  expr_mat.long <- tidyr::gather(expr_mat.long, "sample", "expression", 2:ncol(expr_mat.long))
  # Prevent ggplot auto-ordering 
  expr_mat.long$gene   <- factor(expr_mat.long$gene,   levels = row.names(expr_mat)[hcRow$order])
  expr_mat.long$sample <- factor(expr_mat.long$sample, levels = colnames(expr_mat))
  main.p <- ggplot(expr_mat.long, aes(x=sample, y = gene, fill = expression)) + geom_tile() + theme_minimal_grid()
  main.p <- main.p + main_theme + labs(y = "", x = "") + scale_y_discrete(position = "right", expand = c(0,0)) + theme(axis.text.y = element_text(size = 8), legend.title = element_text(size=7.5), legend.text = element_text(size = 7), legend.position = "right", axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1), axis.ticks.x = element_blank(), plot.margin = unit(c(0,0,0,0.5), "cm"))
  
  if (!showSample) {
    main.p <- main.p + theme(axis.text.x = element_blank())
  }
  if (!showGene) {
    main.p <- main.p + theme(axis.text.y = element_blank())
  }
  if (scale) {
    main.p <- main.p + labs(fill = "Scaled Expression")
    # Rescale and specify where 0 would be 
    rescale_colors <- c("blue","blue2", "blue4", "black", "yellow4", "yellow2", "yellow")
    negative_midpoint <- min(expr_mat.long$expression)/1.5
    positive_midpoint <- max(expr_mat.long$expression)/1.5
    rescale_vals <- scales::rescale(c(min(expr_mat.long$expression), negative_midpoint, -0.5, 0, 0.5, positive_midpoint , max(expr_mat.long$expression)))
    main.p <- main.p + scale_fill_gradientn(colours = rescale_colors, values = rescale_vals)
  } else {
    main.p <- main.p + scale_fill_viridis(option = "B", begin = 0) 
  }
  ##########---------------
  # Plot gene row annot
  # Rearrange annot_df to reflect clustering 
  annot_df <- annot_df[levels(expr_mat.long$gene),]
  main.p <- main.p + theme(axis.text.y = element_text(color = annot_df$Color, size = 6))
  ##########---------------
  # Plot sample col annot
  # Create a subset of samples & Remove tmp column made above
  smple_sub <- smple_df[row.names(smple_df) %in% colnames(expr_mat),]
  smple_sub <- smple_sub[, !(colnames(smple_sub) %in% "tmp"), drop = FALSE]
  smple_df.long <- smple_sub %>% tibble::rownames_to_column(var = "sample")
  smple_df.long <- tidyr::gather(smple_df.long, "groupBy", "Group", 2:ncol(smple_df.long))
  # Prevent ggplot auto-ordering 
  smple_df.long$sample <- factor(smple_df.long$sample, levels =  levels(expr_mat.long$sample))
  ColAnnot.p <- ggplot(smple_df.long, aes(x=sample, y = groupBy, fill = Group)) + geom_tile(color=NA) + scale_fill_manual(values = color_r) + theme_minimal_grid() + main_theme + labs(y="", x="")
  ColAnnot.p <- ColAnnot.p + scale_y_discrete(position = "right", expand = c(0,0)) + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 8), panel.background = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(0.05,0,0,0.5), "cm"), legend.text = element_text(size = 6), legend.title = element_blank())
  ##########---------------
  # ALL TOGETHER NOW
  P <- cowplot::plot_grid(ColAnnot.p, main.p, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.1, 1))
  return(list(actualPlot = P, main = main.p, annot = ColAnnot.p))
} # end do_heatmapSeurat.f
