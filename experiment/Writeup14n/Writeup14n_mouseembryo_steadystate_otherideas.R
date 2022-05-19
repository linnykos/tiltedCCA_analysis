rm(list=ls())
load("../../../out/main/10x_mouseembryo_tiltedcca.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# rna_common <- multiSVD_obj$common_mat_1
# multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
# svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
# multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
# svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
# tmp <- crossprod(svd_2$u, svd_1$u)
# svd_tmp <- svd(tmp)
# rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
# atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
# n <- nrow(rna_common)
# alignment_vec <- sapply(1:n, function(i){
#   df <- data.frame(rna = rna_common[i,],
#                    atac = atac_pred[i,])
#   lm_res <- stats::lm(rna ~ atac, data = df)
#   summary(lm_res)$r.squared
# })

rna_common <- multiSVD_obj$common_mat_1
rna_distinct <- multiSVD_obj$distinct_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.svd_safe(rna_common,
                               check_stability = T,
                               K = ncol(svd_2$u),
                               mean_vec = NULL,
                               rescale = F,
                               scale_max = NULL,
                               sd_vec = NULL)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(rna = rna_common[i,]+rna_distinct[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})

mbrain$alignment <- rank(alignment_vec)
plot1 <-Seurat::FeaturePlot(mbrain, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nAlignment between ATAC and common RNA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_tcca_steadystate-full.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")
