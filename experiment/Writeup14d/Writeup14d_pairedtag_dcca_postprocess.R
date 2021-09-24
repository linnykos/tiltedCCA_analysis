rm(list=ls())
histone_names <- c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3")

cell_types <- c("Astro_Myoc", "Astro_Nnat", "CA1", "CA23", "CGE",
                "CT", "DG", "Endothelial", "Ependymal", "L23",
                "L4", "L5", "L6", "Microglia", "NP",
                "Oligo_MFOL", "Oligo_MOL", "OPC", "PT", "Pvalb",
                "Sst", "Subiculum")
color_vec <- c("#FF6B2C", "#FF8817", "#5E801F", "#06C65B", "#EA6FE6",
               "#294646", "#00803D", "#B26F13", "#804F0F", "#005D7F",
               "#008FC6", "#05A8EB", "#547180", "#EB8E1A", "#487F80",
               "#EB523A", "#B2402D", "#7F2F23", "#6EC6C6", "#803E7E",
               "#B255B0", "#AAEB32")
color_df <- data.frame(celltype = cell_types, color = color_vec)


for(kk in 1:length(histone_names)){
  print(kk)
  
  load(paste0("../../../../out/Writeup14d/Writeup14d_pairedtag_",  
               histone_names[kk], ".RData"))
  
  pairedtag[["common"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[1]], key = "common_", assay = "RNA")
  pairedtag[["distinct"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[2]], key = "distinct_", assay = "RNA")
  pairedtag[["everything"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[3]], key = "everything_", assay = "RNA")
  rm(list = c("rna_embeddings", "rna_frnn")); gc(T)
  
  pairedtag[["common2"]] <- Seurat::CreateDimReducObject(embedding = dna_embeddings[[1]], key = "common2_", assay = "RNA")
  pairedtag[["distinct2"]] <- Seurat::CreateDimReducObject(embedding = dna_embeddings[[2]], key = "distinct2_", assay = "RNA")
  pairedtag[["everything2"]] <- Seurat::CreateDimReducObject(embedding = dna_embeddings[[3]], key = "everything2_", assay = "RNA")
  rm(list = c("dna_embeddings", "dna_frnn")); gc(T)
  
  pairedtag[["combined"]] <- Seurat::CreateDimReducObject(embedding = combined_common_umap, key = "combined_", assay = "RNA")
  rm(list = c("combined_common_umap")); gc(T)
  
  uniq_celltypes <- sort(unique(pairedtag@meta.data$celltype))
  color_vec <- sapply(uniq_celltypes, function(i){
    color_df[which(color_df$celltype == i),"color"]
  })
  
  embedding_1 <- pairedtag[["pca"]]@cell.embeddings
  embedding_2 <- pairedtag[["lsi"]]@cell.embeddings
  
  svd_1 <- irlba::irlba(embedding_1, nv = 3)
  embedding_1 <- embedding_1/svd_1$d[1]*nrow(embedding_1)
  svd_2 <- irlba::irlba(embedding_2, nv = 3)
  embedding_2 <- embedding_2/svd_2$d[1]*nrow(embedding_2)
  
  set.seed(10)
  both_embedding <- Seurat::RunUMAP(cbind(embedding_1, embedding_2))
  both_embedding <- both_embedding@cell.embeddings
  pairedtag[["both"]] <- Seurat::CreateDimReducObject(embedding = both_embedding, key = "both_", assay = "RNA")
  
  ###########################
  
  # done with all the calculations (for now). now plot
  
  title_vec <- c("Common view", "Distinct view", "Everything view")
  main_vec <- c("common", "distinct", "everything")
  
  # plot RNA embeddings
  for(i in 1:3){
    plot1 <- Seurat::DimPlot(pairedtag, reduction = main_vec[i],
                             group.by = "celltype", label = TRUE,
                             cols = color_vec,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ", histone_names[kk], " (RNA)\n", title_vec[i]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_pairedtag_",  histone_names[kk], "_dcca_rna_", main_vec[i], "_umap.png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
  
  # plot DNA embeddings
  for(i in 1:3){
    plot1 <- Seurat::DimPlot(pairedtag, reduction = paste0(main_vec[i], "2"),
                             group.by = "celltype", label = TRUE,
                             cols = color_vec,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ", histone_names[kk], " (Histone)\n", title_vec[i]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_pairedtag_",  histone_names[kk], "_dcca_dna_", main_vec[i], "_umap.png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
  
  # plot the combined view
  plot1 <- Seurat::DimPlot(pairedtag, reduction = "combined",
                           group.by = "celltype", label = TRUE,
                           cols = color_vec,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ", histone_names[kk], "\nBoth Common"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_pairedtag_",  histone_names[kk], "_dcca_both_common_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  # plot the everything combined for both
  plot1 <- Seurat::DimPlot(pairedtag, reduction = "both",
                           group.by = "celltype", label = TRUE,
                           cols = color_vec,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ", histone_names[kk], "\nBoth Everything"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_pairedtag_",  histone_names[kk], "_dcca_both_everything_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

################################

for(kk in 1:length(histone_names)){
  print(kk)
  
  load(paste0("../../../../out/Writeup14d/Writeup14d_pairedtag_",  
              histone_names[kk], ".RData"))
  
  png(paste0("../../../../out/figures/Writeup14d/Writeup14d_pairedtag_",
             histone_names[kk], "_dcca_summary.png"),
      height = 1500, width = 1500, units = "px", res = 300)
  par(mar = c(5,5,5,5))
  multiomicCCA::plot_summary(dcca_res,
                             main = paste0("Paired-Tag: ", histone_names[kk]))
  graphics.off()
  
  png(paste0("../../../../out/figures/Writeup14d/Writeup14d_pairedtag_",
      histone_names[kk], "_dcca_scores.png"),
      height = 1200, width = 2500, units = "px", res = 300)
  par(mfrow = c(1,3), mar = c(4,4,4,0.5))
  multiomicCCA::plot_scores_heatmap.dcca(dcca_res,
                                         membership_vec = as.factor(pairedtag@meta.data$celltype),
                                         log_scale = T, scaling_power = 2)
  graphics.off()
}