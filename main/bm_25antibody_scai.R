rm(list=ls())
library(Seurat)

load("../../../out/main/citeseq_bm25_preprocessed.RData")
source("bm_25antibody_colorPalette.R")
source("scai.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

####

Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@data[Seurat::VariableFeatures(object = bm),])
Seurat::DefaultAssay(bm) <- "ADT"
mat_2 <- Matrix::t(bm[["ADT"]]@data)

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

mat_1b <- as.matrix(mat_1b)
mat_2b <- as.matrix(mat_2b)

set.seed(10)
scai_res <- scai(mat_1 = mat_1b, 
                 mat_2 = mat_2b, 
                 r = 30, 
                 gamma = 0, 
                 max_iter = 100)

save(scai_res, bm,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_bm25_scAI.RData")
