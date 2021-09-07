make_all_embedding_plots <- function(seurat_obj, 
                                     title_all,
                                     title_vec, 
                                     main_vec, 
                                     file_prefix,
                                     file_suffix){
  # plot RNA embedding
  for(i in 1:3){
    plot1 <- Seurat::DimPlot(seurat_obj, reduction = main_vec[i], 
                             group.by = "celltype", label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0(title_all, " (RNA)\n", title_vec[i]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0(file_prefix, "_rna_", main_vec[i], "_umap", file_suffix, ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
  
  # plot ATAC embeddings
  for(i in 1:3){
    plot1 <- Seurat::DimPlot(seurat_obj, reduction = paste0(main_vec[i], "2"), 
                             group.by = "celltype", label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0(title_all, " (ATAC)\n", title_vec[i])) 
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0(file_prefix, "_atac_", main_vec[i], "_umap", file_suffix, ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
  
  # plot the combined view
  plot1 <- Seurat::DimPlot(seurat_obj, reduction = "combined", 
                           group.by = "celltype", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(title_all, "\nBoth Common") 
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0(file_prefix, "_common_umap", file_suffix, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  # plot the everything combined for both
  plot1 <- Seurat::DimPlot(seurat_obj, reduction = "both", 
                           group.by = "celltype", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(title_all, "\nBoth Everything") 
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0(file_prefix, "_both_everything_umap", file_suffix, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")

  invisible()
}