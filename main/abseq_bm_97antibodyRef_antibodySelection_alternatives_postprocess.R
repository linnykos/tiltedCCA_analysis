rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/abseq_bm97Ref_varSelect_alternatives.RData")
source("bm_97antibodyRef_colorPalette.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

for(i in 1:3){
  print(i)
  reduction_name <- paste0("adt.umap", i)
  plot1 <- Seurat::DimPlot(bm, reduction = reduction_name,
                           group.by = "ct", label = TRUE,
                           repel = TRUE, label.size = 2.5,
                           cols = col_palette,
                           raster = FALSE)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nADT: Selected antibodies: ", names(panel_list)[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_varSelect-alternative", i, "_adt-umap.png"),
                  plot1, device = "png", width = 11, height = 5, units = "in")
  
  reduction_name <- paste0("consensusUMAP", i)
  plot1 <- Seurat::DimPlot(bm, reduction = reduction_name,
                           group.by = "ct", label = TRUE,
                           repel = TRUE, label.size = 2.5,
                           cols = col_palette,
                           raster = FALSE)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nConsensus PCA: Selected antibodies: ", names(panel_list)[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_varSelect-alternative", i, "_consensusPCA-umap.png"),
                  plot1, device = "png", width = 11, height = 5, units = "in")
}


for(i in 1:4){
  print(i)
  reduction_name <- paste0("adt.umap", i)
  plot1 <- Seurat::DimPlot(bm, reduction = reduction_name,
                           group.by = "ct",
                           cols = col_palette,
                           raster = FALSE)
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_varSelect-alternative", i, "_adt-umap_cleaned.png"),
                  plot1, device = "png", width = 3, height = 3, units = "in",
                  dpi = 500)
  
  
  reduction_name <- paste0("consensusUMAP", i)
  plot1 <- Seurat::DimPlot(bm, reduction = reduction_name,
                           group.by = "ct",
                           cols = col_palette,
                           raster = FALSE)
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_varSelect-alternative", i, "_consensusPCA-umap_cleaned.png"),
                  plot1, device = "png", width = 3, height = 3, units = "in",
                  dpi = 500)
}


