rm(list=ls())
library(Seurat)
library(tiltedCCA)

load("../../../out/main/citeseq_bm25_preprocessed.RData")
source("bm_25antibody_colorPalette.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# consensus pca
Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@data[Seurat::VariableFeatures(object = bm),])
Seurat::DefaultAssay(bm) <- "ADT"
mat_2 <- Matrix::t(bm[["ADT"]]@data[Seurat::VariableFeatures(object = bm),])

mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

##################

dims_1 = 1:30
dims_2 = 1:18
dims_consensus = 1:30
center_1 = T
center_2 = T
recenter_1 = F
recenter_2 = F
rescale_1 = F
rescale_2 = F
scale_1 = T
scale_2 = T
verbose = 1
scale_max_1 = NULL
scale_max_2 = NULL
normalize_row = F
normalize_singular_value = T
center_consensus = T
recenter_consensus = F
rescale_consensus = F
scale_consensus = T
scale_max_consensus = NULL

svd_1 <- tiltedCCA:::.get_SVD(center = center_1, input_obj = mat_1,
                              dims = dims_1, scale = scale_1, scale_max = scale_max_1)
svd_2 <- tiltedCCA:::.get_SVD(center = center_2, input_obj = mat_2,
                              dims = dims_2, scale = scale_2, scale_max = scale_max_2)

dimred_1 <- tiltedCCA:::.normalize_svd(input_obj = svd_1,
                                       averaging_mat = NULL,
                                       normalize_row = normalize_row,
                                       normalize_singular_value = normalize_singular_value,
                                       recenter = recenter_1,
                                       rescale = rescale_1)
dimred_2 <- tiltedCCA:::.normalize_svd(input_obj = svd_2,
                                       averaging_mat = NULL,
                                       normalize_row = normalize_row,
                                       normalize_singular_value = normalize_singular_value,
                                       recenter = recenter_2,
                                       rescale = rescale_2,
                                       tol = 1e-4)

if(verbose > 0) print("Computing Consensus PCA")
stopifnot(nrow(dimred_1) == nrow(dimred_2))
dimred_combined <- cbind(dimred_1, dimred_2)

svd_consensus <- tiltedCCA:::.get_SVD(center = center_consensus, 
                                      input_obj = dimred_combined,
                                      dims = dims_consensus, 
                                      scale = scale_consensus, 
                                      scale_max = scale_max_consensus)
dim(svd_consensus$v)
rna_loadings <- apply(svd_consensus$v[1:30,1:18], 2, function(x){mean(abs(x))})
adt_loadings <- apply(svd_consensus$v[41:48,1:18], 2, function(x){mean(abs(x))})

cbind(round(100*rna_loadings,2), round(100*adt_loadings,2))
