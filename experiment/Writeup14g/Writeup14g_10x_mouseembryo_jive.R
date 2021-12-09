rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")
source("../../simulation/jive.R")
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)
jive_res <- jive(mat_1, mat_2, r = 30)
save.image(file = "../../../../out/Writeup14g/Writeup14g_citeseq_bm25_jive.RData")