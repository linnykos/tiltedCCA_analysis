rm(list=ls())
library(Seurat)
library(Signac)

histone_names <- c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3")

# create the color vector
cell_types <- sort(unique(unlist(lapply(1:length(histone_names), function(i){
  load(paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/seurat_",
              histone_names[i], ".RData"))
  unique(pairedtag@meta.data$celltype)
}))))
# [1] "Astro_Myoc"  "Astro_Nnat"  "CA1"         "CA23"        "CGE"        
# [6] "CT"          "DG"          "Endothelial" "Ependymal"   "L23"        
# [11] "L4"          "L5"          "L6"          "Microglia"   "NP"         
# [16] "Oligo_MFOL"  "Oligo_MOL"   "OPC"         "PT"          "Pvalb"      
# [21] "Sst"         "Subiculum"
color_vec <- c("#FF6B2C", "#FF8817", "#5E801F", "#06C65B", "#EA6FE6",
               "#294646", "#00803D", "#B26F13", "#804F0F", "#005D7F",
               "#008FC6", "#05A8EB", "#547180", "#EB8E1A", "#487F80",
               "#EB523A", "#B2402D", "#7F2F23", "#6EC6C6", "#803E7E",
               "#B255B0", "#AAEB32")
color_df <- data.frame(celltype = cell_types, color = color_vec)

for(i in 1:length(histone_names)){
  print(i)
  
  load(paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/seurat_",
              histone_names[i], ".RData"))
  print(quantile(pairedtag[["DNA"]]@counts@x))
  
  uniq_celltypes <- sort(unique(pairedtag@meta.data$celltype))
  color_vec <- sapply(uniq_celltypes, function(i){
    color_df[which(color_df$celltype == i),"color"]
  })
  
  plot1 <- Seurat::DimPlot(pairedtag, reduction = "umap",
                           group.by = "celltype", label = TRUE,
                           cols = color_vec,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ",  histone_names[i]," (RNA)"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_pairedtag_",
                                    histone_names[i],
                                    "_rna.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  plot1 <- Seurat::DimPlot(pairedtag, reduction = "umap.dna",
                           group.by = "celltype", label = TRUE,
                           cols = color_vec,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ",  histone_names[i]," (Histone)"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_pairedtag_",
                                    histone_names[i],
                                    "_histone.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  plot1 <- Seurat::DimPlot(pairedtag, reduction = "wnn.umap",
                           group.by = "celltype", label = TRUE,
                           cols = color_vec,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ",  histone_names[i]," (WNN)"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_pairedtag_",
                                    histone_names[i],
                                    "_wnn.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}
