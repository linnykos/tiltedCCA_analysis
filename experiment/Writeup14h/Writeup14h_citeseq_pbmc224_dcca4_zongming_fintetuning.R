rm(list=ls())
load("../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_dcca4_zongming.RData")
library(Seurat); library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

dcca2 <- multiomicCCA:::fine_tuning(dcca_res, 
                                   max_iter = 5,
                                   fix_tilt_perc = NA,
                                   temp_path = "../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_dcca4_zongming_tmp.RData",
                                   verbose = T)

save(pbmc, dcca_res, dcca_res2, 
     rank_1, rank_2, nn, date_of_run, session_info,
     file = "../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_dcca4_zongming_finetuning.RData")


