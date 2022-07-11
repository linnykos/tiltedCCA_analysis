rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_mouseembryo_tiltedcca.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

cell_idx <- unique(c(which(mbrain$Lineage1 == 1), which(mbrain$Lineage2 == 1)))
input_obj = multiSVD_obj
bool_use_denoised = F
bool_include_intercept = T
bool_use_metacells = F
cor_threshold = 0.8
input_assay = 1
num_variables = 50
sd_quantile = 0.75
seurat_obj = mbrain
seurat_assay = "SCT"
seurat_slot = "data"

res <- tiltedCCA:::postprocess_graph_alignment(
  input_obj = input_obj,
  bool_use_denoised = bool_use_denoised,
  bool_include_intercept = bool_include_intercept,
  bool_use_metacells = bool_use_metacells,
  cell_idx = cell_idx,
  input_assay = input_assay,
  return_everything_mat = T,
  seurat_obj = seurat_obj,
  seurat_assay = seurat_assay,
  seurat_slot = seurat_slot,
  verbose = 2
)

png("../../../out/figures/main/10x_mbrain_tcca_develompentalGenes_histogram.png",
    height = 1100*.5, width = 950*.5, res = 500, units = "px")
par(mar = c(1,1,0,0))
hist(res$alignment, col = rgb(55, 100, 234, maxColorValue = 255),
     xlab = "", ylab = "", main = "")
graphics.off()