make_all_embedding_plots <- function(seurat_obj, 
                                     title_all,
                                     title_vec, 
                                     main_vec, 
                                     file_prefix,
                                     file_suffix,
                                     group.by = "celltype"){
  # plot RNA embedding
  for(i in 1:3){
    plot1 <- Seurat::DimPlot(seurat_obj, reduction = main_vec[i], 
                             group.by = group.by, label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0(title_all, " (RNA)\n", title_vec[i]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0(file_prefix, "_rna_", main_vec[i], "_umap", file_suffix, ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
  
  # plot ATAC embeddings
  for(i in 1:3){
    plot1 <- Seurat::DimPlot(seurat_obj, reduction = paste0(main_vec[i], "2"), 
                             group.by = group.by, label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0(title_all, " (ATAC)\n", title_vec[i])) 
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0(file_prefix, "_atac_", main_vec[i], "_umap", file_suffix, ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
  
  # plot the combined view
  plot1 <- Seurat::DimPlot(seurat_obj, reduction = "combined", 
                           group.by = group.by, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(title_all, "\nBoth Common") 
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0(file_prefix, "_common_umap", file_suffix, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  # plot the everything combined for both
  plot1 <- Seurat::DimPlot(seurat_obj, reduction = "both", 
                           group.by = group.by, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(title_all, "\nBoth Everything") 
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0(file_prefix, "_both_everything_umap", file_suffix, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")

  invisible()
}



plot_graph <- function(embedding, 
                       clustering,
                       adj_mat,
                       feature_vec = NA,
                       zlim = range(feature_vec),
                       base_color = grDevices::rgb(0.49, 0, 0.81),
                       missing_color = "gray70",
                       col_vec = scales::hue_pal()(nrow(adj_mat)),
                       asp = T,
                       bins = 100,
                       cex_missing = 0.5,
                       cex_normal = 1,
                       ...){
  selected_clusters <- rownames(adj_mat)
  stopifnot(ncol(embedding) == 2,
            length(col_vec) == length(selected_clusters))
  if(!all(is.na(feature_vec))) stopifnot(length(feature_vec) == length(selected_clusters))
  
  mat <- .form_metacell_matrix(embedding, clustering)
  
  if(all(!is.na(feature_vec))){
    stopifnot(length(feature_vec) == nrow(adj_mat))
    col_ramp <- grDevices::colorRampPalette(c("white", base_color))(bins)
    val_vec <- seq(zlim[1], zlim[2], length = bins)
    col_idx <- sapply(feature_vec, function(x){
      which.min(abs(x - val_vec))
    })
    col_vec <- col_ramp[col_idx]
  } 
  
  n <- nrow(mat)
  full_col_vec <- rep(missing_color, n)
  for(i in 1:length(selected_clusters)){
    full_col_vec[which(rownames(mat) == selected_clusters[i])] <- col_vec[i]
  }
  cluster_mapping <- sapply(selected_clusters, function(x){
    which(rownames(mat) == x)
  })
 
  graphics::plot(NA, 
                 xlim = range(mat[,1]),
                 ylim = range(mat[,2]),
                 asp = asp,
                 ...)
  
  if(any(missing_color %in% full_col_vec)){
    idx <- which(full_col_vec == missing_color)
    
    for(k in idx){
      points(mat[k,1], mat[k,2], 
             pch = 16, 
             col = full_col_vec[k],
             cex = cex_missing)
    }
    
  } 
  
  for(j in 1:nrow(adj_mat)){
    idx <- which(adj_mat[j,] != 0)
    for(j2 in idx){
      lines(mat[cluster_mapping[c(j,j2)],])
    }
  }
  
  for(k in cluster_mapping){
    points(mat[k,1], mat[k,2], 
           pch = 16, 
           col = "black",
           cex = cex_normal+1)
    points(mat[k,1], mat[k,2], 
           pch = 16, 
           col = full_col_vec[k],
           cex = cex_normal)
  }
  
  invisible()
}