rm(list=ls())
load("../../../../out/Writeup14k/Writeup14k_citeseq_bm25_preprocessed.RData")
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca2.RData")
load("../../../../out/Writeup14i/Writeup14i_citeseq_pbmc224_dcca_variable_full_de.RData")
library(Seurat)
source("variable_distinct_selection.R")

adt_mat <- tcrossprod(multiomicCCA:::.mult_mat_vec(dcca_res2$svd_2$u, dcca_res2$svd_2$d),
                      dcca_res2$svd_2$v)

adt_correlation <- sapply(1:ncol(adt_mat), function(j){
  df <- data.frame(dcca_res2$svd_1$u)
  df <- cbind(df, adt_mat[,j])
  colnames(df)[ncol(df)] <- "y"
  
  lm_res <- stats::lm(y ~ ., data = df)
  summary(lm_res)$r.squared
})

adt_correlation2 <- sapply(1:ncol(adt_mat), function(j){
  df <- data.frame(dcca_res2$common_score)
  df <- cbind(df, adt_mat[,j])
  colnames(df)[ncol(df)] <- "y"
  
  lm_res <- stats::lm(y ~ ., data = df)
  summary(lm_res)$r.squared
})

###############################

bm <- SeuratData::LoadData(ds = "bmcite")
cd_idx <- grep("^CD", rownames(bm[["RNA"]]))
sort(rownames(bm[["RNA"]])[cd_idx])

