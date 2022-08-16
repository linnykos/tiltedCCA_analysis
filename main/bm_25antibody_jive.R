rm(list=ls())
library(Seurat)
library(r.jive)

load("../../../out/main/citeseq_bm25_preprocessed.RData")
source("bm_25antibody_colorPalette.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

####

Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@scale.data[Seurat::VariableFeatures(object = bm),])
Seurat::DefaultAssay(bm) <- "ADT"
mat_2 <- Matrix::t(bm[["ADT"]]@scale.data)

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

data_list <- list(RNA = t(mat_1b), ADT = t(mat_2b))

set.seed(10)
jive_results <- r.jive::jive(data = data_list,
                             rankJ = 30,
                             rankA = c(30, 18),
                             method = "given",
                             showProgress = TRUE)

save(jive_results, bm,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_bm25_JIVE.RData")