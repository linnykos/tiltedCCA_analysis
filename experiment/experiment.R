rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

cell_idx <- unique(c(which(!is.na(greenleaf$Lineage1)), which(!is.na(greenleaf$Lineage2))))

input_obj = multiSVD_obj
bool_use_denoised = F
bool_include_intercept = T
bool_use_metacells = F
cell_idx = cell_idx
cor_threshold = 0.8
input_assay = 1
num_variables = 50
seurat_obj = greenleaf
seurat_assay = "SCT"
seurat_slot = "data"
verbose = 2

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
  verbose = verbose
)

alignment_vec <- res$alignment
everything_mat <- res$everything_mat
if(any(is.na(alignment_vec))){
  everything_mat <- everything_mat[,which(!is.na(alignment_vec)), drop = F]
  alignment_vec <- alignment_vec[which(!is.na(alignment_vec))]
}
stopifnot(all(names(alignment_vec) == names(everything_mat)))
sd_vec <- matrixStats::colSds(everything_mat)
p <- length(alignment_vec)

