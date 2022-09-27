rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_mouseembryo_tiltedcca_RNA-geneActivity.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

cell_idx <- unique(c(which(mbrain$Lineage1 == 1), which(mbrain$Lineage2 == 1)))
input_obj = multiSVD_obj
bool_use_denoised = T
bool_include_intercept = T
bool_use_metacells = F
cor_threshold = 0.95
num_variables = 50
sd_quantile = 0.05
seurat_obj = mbrain
seurat_slot = "data"

res <- tiltedCCA:::postprocess_graph_alignment(
  input_obj = input_obj,
  bool_use_denoised = bool_use_denoised,
  bool_include_intercept = bool_include_intercept,
  bool_use_metacells = bool_use_metacells,
  cell_idx = cell_idx,
  input_assay = 1,
  return_everything_mat = T,
  seurat_obj = seurat_obj,
  seurat_assay = "SCT",
  seurat_slot = seurat_slot,
  verbose = 2
)

res2 <- tiltedCCA:::postprocess_graph_alignment(
  input_obj = input_obj,
  bool_use_denoised = bool_use_denoised,
  bool_include_intercept = bool_include_intercept,
  bool_use_metacells = bool_use_metacells,
  cell_idx = cell_idx,
  input_assay = 2,
  return_everything_mat = T,
  seurat_obj = seurat_obj,
  seurat_assay = "geneActivity",
  seurat_slot = seurat_slot,
  verbose = 2
)

xlim <- range(c(res$alignment[!is.na(res$alignment)], res2$alignment[!is.na(res2$alignment)]))
breaks <- seq(xlim[1]-.05, xlim[2]+.05, length.out = 15)

png("../../../out/figures/main/10x_mbrain_tcca_develompentalGenes_histogram.png",
    height = 850*.5, width = 950*.5, res = 500, units = "px")
par(mar = c(1,1,0.5,0.5))
hist(res$alignment[!is.na(res$alignment)], breaks = breaks,
     col = rgb(55, 100, 234, maxColorValue = 255),
     xlab = "", ylab = "", main = "", ylim = c(0, 500))
graphics.off()


png("../../../out/figures/main/10x_mbrain_tcca_develompentalGenes_histogram2.png",
    height = 850*.5, width = 950*.5, res = 500, units = "px")
par(mar = c(1,1,0.5,0.5))
hist(res2$alignment[!is.na(res2$alignment)], 
     breaks = breaks, col = rgb(191, 74, 223, maxColorValue = 255),
     xlab = "", ylab = "", main = "",  ylim = c(0, 500))
graphics.off()