rm(list=ls())
library(Seurat)
library(Signac)

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

cell_list <- list(Non_neurons = c("Astro_Myoc", "Astro_Nnat", "OPC", "Microglia", 
                                  "Oligo_MFOL", "Oligo_MOL", "Endothelial", "Ependymal"),
                  Inhibitory_neurons = c("Pvalb", "Sst", "CGE"),
                  Cortical_neurons = c("CT", "NP", "L6", "L5", "L4", "L23", "PT"),
                  Hippocampal_neurons = c("DG", "Subiculum", "CA1", "CA23"))

for(i in 1:length(histone_names)){
  print(i)
  
  load(paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/seurat_",
              histone_names[i], ".RData"))
  print(quantile(pairedtag[["DNA"]]@counts@x))
  cellgroup <- rep(NA, nrow(pairedtag@meta.data))
  for(kk in 1:length(cell_list)){
    cellgroup[which(pairedtag@meta.data[,"celltype"] %in% cell_list[[kk]])] <- names(cell_list)[kk]
  }
  pairedtag[["cellgroup"]] <- cellgroup
  
  uniq_celltypes <- sort(unique(pairedtag@meta.data$celltype))
  color_vec <- sapply(uniq_celltypes, function(i){
    color_df[which(color_df$celltype == i),"color"]
  })
  
  # plot umap
  plot1 <- Seurat::DimPlot(pairedtag, reduction = "umap.dna",
                           group.by = "celltype", label = TRUE,
                           cols = color_vec,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ",  histone_names[i]," (Histone)"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_pairedtag_",
                                    histone_names[i],
                                    "_histone.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  # plot the library size
  plot1 <- Seurat::FeaturePlot(pairedtag, reduction = "umap.dna",
                               features = "nCount_DNA")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ",  histone_names[i]," (Histone)\nnCount_DNA"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_pairedtag_",
                                    histone_names[i],
                                    "_histone_nCount_DNA.png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
  # plot each cell group separately
  
  p_list <- lapply(sort(unique(pairedtag@meta.data[,"cellgroup"])), function(cellgroup){
    p0 <- Seurat::DimPlot(pairedtag, cells.highlight = rownames(pairedtag@meta.data)[which(pairedtag@meta.data[,"cellgroup"] == cellgroup)])
    p0 <- p0 + ggplot2::theme(legend.position="none") + ggplot2::ggtitle(paste0("Group: ", cellgroup))
    p0
  })
  
  p <- cowplot::plot_grid(p_list[[1]], p_list[[2]], p_list[[3]], p_list[[4]])
  cowplot::save_plot(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_pairedtag_",
                                       histone_names[i],
                                       "_histone_separate.png"),
                     p, ncol = 2, nrow = 2, base_asp = 1, device = "png")
  
  
}
