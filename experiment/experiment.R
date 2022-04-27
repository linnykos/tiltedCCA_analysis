rm(list=ls())
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")

library(Seurat); library(Signac)

common_mat <- tiltedCCA:::.convert_common_score_to_mat(common_score = multiSVD_obj$tcca_obj$common_score,
                                                       score_1 = multiSVD_obj$cca_obj$score_1,
                                                       score_2 = multiSVD_obj$cca_obj$score_2,
                                                       svd_1 = multiSVD_obj$svd_1,
                                                       svd_2 = multiSVD_obj$svd_2)

set.seed(10)
seurat_umap <- Seurat::RunUMAP(common_mat, 
                               assay = "SCT")
rownames(seurat_umap@cell.embeddings) <- rownames(score_1)

pbmc[["tmp"]] <- Seurat::CreateDimReducObject(seurat_umap@cell.embeddings,
                                              assay = "SCT")

plot1 <- Seurat::DimPlot(pbmc, reduction = "tmp",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/Writeup14n/citeseq_pbmc224_tcca-umap_common_version1.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

############################################

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = F,
                                                    bool_modality_2_full = F)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
dimred_1 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, what = "common_dimred")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
dimred_2 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, what = "common_dimred")
dimred <- cbind(dimred_1, dimred_2)

set.seed(10)
seurat_umap <- Seurat::RunUMAP(dimred, 
                               assay = "SCT")
rownames(seurat_umap@cell.embeddings) <- rownames(score_1)

pbmc[["tmp"]] <- Seurat::CreateDimReducObject(seurat_umap@cell.embeddings,
                                              assay = "SCT")

plot1 <- Seurat::DimPlot(pbmc, reduction = "tmp",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/Writeup14n/citeseq_pbmc224_tcca-umap_common_version2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

#####################################

multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
dimred_1b <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, what = "common_mat")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
dimred_2b <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, what = "common_mat")
dimredb <- cbind(dimred_1b, dimred_2b)



set.seed(10)
seurat_umap <- Seurat::RunUMAP(dimredb, 
                               assay = "SCT")
rownames(seurat_umap@cell.embeddings) <- rownames(score_1)

pbmc[["tmp"]] <- Seurat::CreateDimReducObject(seurat_umap@cell.embeddings,
                                              assay = "SCT")
plot1 <- Seurat::DimPlot(pbmc, reduction = "tmp",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/Writeup14n/citeseq_pbmc224_tcca-umap_common_version3.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##################

tmp <- multiSVD_obj$common_mat_1
param <- tiltedCCA:::.get_param(multiSVD_obj)
normalize_row <- param$svd_normalize_row
normalize_singular_value <- param$svd_normalize_singular_value
recenter <- param$svd_recenter_1
rescale <- param$svd_rescale_1
center <- param$svd_center_1
dims <- param$svd_dims_1; dims <- dims - min(dims) + 1
dims <- pmin(dims, ncol(tmp))
scale <- param$svd_scale_1

dimred_tmp <- tiltedCCA:::.get_Dimred(input_obj = multiSVD_obj, 
                                      normalize_singular_value = normalize_singular_value,
                                      center = center,
                                      dims = dims,
                                      scale = scale)

if(normalize_row){
  l2_vec <- apply(dimred_tmp, 1, function(x){tiltedCCA:::.l2norm(x)})
  l2_vec[l2_vec <= 1e-4] <- 1e-4
  dimred_tmp <- tiltedCCA:::.mult_vec_mat(1/l2_vec, dimred_tmp)
}

###

multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
zz <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "common_mat")
yy <- tiltedCCA:::.get_Dimred(input_obj = zz, 
                              normalize_singular_value = normalize_singular_value,
                              center = center,
                              dims = dims,
                              scale = scale) # this is the problematic part...

yy2 <- tiltedCCA:::.get_Dimred(input_obj = zz, 
                              normalize_singular_value = normalize_singular_value,
                              center = F,
                              dims = dims,
                              scale = F)
yy3b <- tiltedCCA:::.get_SVD(input_obj = multiSVD_obj, 
                            center = center,
                            dims = dims,
                            scale = scale)
n <- nrow(yy3b$u)
if(normalize_singular_value) yy3b$d <- yy3b$d*sqrt(n)/yy3b$d[1]
yy4b <- tiltedCCA:::.mult_mat_vec(yy3b$u, yy3b$d)

##

yy3 <- tiltedCCA:::.get_SVD(input_obj = zz, 
                            center = center,
                            dims = dims,
                            scale = scale)
n <- nrow(yy3$u)
if(normalize_singular_value) yy3$d <- yy3$d*sqrt(n)/yy3$d[1]
yy4 <- tiltedCCA:::.mult_mat_vec(yy3$u, yy3$d)

