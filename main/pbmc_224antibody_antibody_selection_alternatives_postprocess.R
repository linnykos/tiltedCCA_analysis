rm(list=ls())
load("../../../out/main/citeseq_pbmc224_varSelect_alternatives.RData")
source("pbmc_224antibody_colorPalette.R")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

for(i in 1:3){
  print(i)
  reduction_name <- paste0("adt.umap", i)
  plot1 <- Seurat::DimPlot(pbmc, reduction = reduction_name,
                           group.by = "celltype.l2", label = TRUE,
                           repel = TRUE, label.size = 2.5,
                           cols = col_palette,
                           raster = FALSE)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nADT: Selected antibodies: ", names(panel_list)[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_varSelect-alternative", i, "_adt-umap_cleaned.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  reduction_name <- paste0("consensusPCA", i)
  plot1 <- Seurat::DimPlot(pbmc, reduction = reduction_name,
                           group.by = "celltype.l2", label = TRUE,
                           repel = TRUE, label.size = 2.5,
                           cols = col_palette,
                           raster = FALSE)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nConsensus PCA: Selected antibodies: ", names(panel_list)[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_varSelect-alternative", i, "_consensusPCA-umap_cleaned.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

