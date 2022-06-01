rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

cell_idx <- unique(c(which(!is.na(greenleaf$Lineage1)), which(!is.na(greenleaf$Lineage2))))
selection_res <- tiltedCCA:::postprocess_smooth_variable_selection(
  input_obj = multiSVD_obj,
  bool_use_denoised = F,
  bool_include_intercept = T,
  bool_use_metacells = F,
  cell_idx = cell_idx,
  cor_threshold = 0.8,
  input_assay = 1,
  num_variables = 50,
  seurat_obj = greenleaf,
  seurat_assay = "SCT",
  seurat_slot = "data",
  verbose = 2
)

Seurat::DefaultAssay(greenleaf) <- "SCT"
plot1 <- Seurat::FeaturePlot(greenleaf, 
                             features = selection_res$selected_variables[1:25],
                             reduction = "common_tcca",
                             ncol = 5)
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_develompentalGenes1.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")

plot1 <- Seurat::FeaturePlot(greenleaf, 
                             features = selection_res$selected_variables[26:50],
                             reduction = "common_tcca",
                             ncol = 5)
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_develompentalGenes2.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")

save(selection_res, date_of_run, session_info, cell_idx,
     file = "../../../out/main/10x_greenleaf_developmentalGenes.RData")