rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)
source("bm_25antibody_colorPalette.R")

percentage_vec <- c(0.25, 0.5, 1) # c(0.1, 0.25, 0.5, 1)
for(percentage in percentage_vec){
  print(paste0("Trying percentage ", percentage))
  load(paste0("../../../out/main/citeseq_bm25_tcca_added-noise_", percentage, ".RData"))
  
  print("Plotting")
  plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap",
                           group.by = "celltype.l2", 
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_added-noise_", percentage, "_rna-umap_cleaned.png"),
                  plot1, device = "png", width = 3, height = 3, units = "in",
                  dpi = 500)
  
  plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                           group.by = "celltype.l2", 
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_added-noise_", percentage, "_adt-umap_cleaned.png"),
                  plot1, device = "png", width = 3, height = 3, units = "in",
                  dpi = 500)
  
  plot1 <- Seurat::DimPlot(bm, reduction = "common_tcca",
                           group.by = "celltype.l2", 
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_added-noise_", percentage, "_tcca-umap_common_cleaned.png"),
                  plot1, device = "png", width = 3, height = 3, units = "in",
                  dpi = 500)
  
  plot1 <- Seurat::DimPlot(bm, reduction = "distinct1_tcca",
                           group.by = "celltype.l2", 
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_added-noise_", percentage, "_tcca-umap_distinct1_cleaned.png"),
                  plot1, device = "png", width = 3, height = 3, units = "in",
                  dpi = 500)
  
  plot1 <- Seurat::DimPlot(bm, reduction = "distinct2_tcca",
                           group.by = "celltype.l2", 
                           cols = col_palette)
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_added-noise_", percentage, "_tcca-umap_distinct2_cleaned.png"),
                  plot1, device = "png", width = 3, height = 3, units = "in",
                  dpi = 500)
}
