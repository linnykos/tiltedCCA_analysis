rm(list=ls())
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")

library(Seurat); library(Signac)

input_obj <- multiSVD_obj
seurat_assay <- "SCT"

input_obj$param$svd_center_1 <- F
input_obj$param$svd_center_2 <- T
input_obj$param$svd_scale_1 <- F
input_obj$param$svd_scale_2 <- F

input_obj <- tiltedCCA:::.set_defaultAssay(input_obj, assay = 2)
dimred <- tiltedCCA:::.get_tCCAobj(input_obj, apply_postDimred = T, what = "distinct_mat")
seurat_umap <- Seurat::RunUMAP(dimred, 
                               assay = seurat_assay)
rownames(seurat_umap@cell.embeddings) <- rownames(dimred)
pbmc[["tmp"]] <- Seurat::CreateDimReducObject(seurat_umap@cell.embeddings,
                                              assay = seurat_assay)

plot1 <- Seurat::DimPlot(pbmc, reduction = "tmp",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/citeseq_pbmc224_tcca-umap_distinct_version2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##########################

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = F,
                                                    bool_modality_2_full = F)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
dimred <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "distinct_dimred")

seurat_assay <- "SCT"
seurat_umap <- Seurat::RunUMAP(dimred, 
                               assay = seurat_assay)
rownames(seurat_umap@cell.embeddings) <- rownames(dimred)
pbmc[["tmp"]] <- Seurat::CreateDimReducObject(seurat_umap@cell.embeddings,
                                              assay = seurat_assay)

plot1 <- Seurat::DimPlot(pbmc, reduction = "tmp",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/citeseq_pbmc224_tcca-umap_distinct_version3.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

############

tmp <- multiSVD_obj$common_mat_1
mean_vec <- Matrix::colMeans(tmp)
sd_vec <- matrixStats::colSds(tmp)
irlba_tmp <- irlba::irlba(tmp, nv = 5, center = mean_vec, scale = sd_vec)
irlba_tmp$d

param <- tiltedCCA:::.get_param(multiSVD_obj)
normalize_row <- param$svd_normalize_row
normalize_singular_value <- param$svd_normalize_singular_value
recenter <- param$svd_recenter_1
rescale <- param$svd_rescale_1
center <- param$svd_center_1
dims <- param$svd_dims_1; dims <- dims - min(dims) + 1
dims <- pmin(dims, ncol(tmp))
scale <- param$svd_scale_1

multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
zz <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "common_mat")
yy <- tiltedCCA:::.get_Dimred(input_obj = zz, 
                              normalize_singular_value = normalize_singular_value,
                              center = center,
                              dims = dims,
                              scale = scale) # this is the problematic part...
yy3 <- tiltedCCA:::.get_SVD(input_obj = zz, 
                            center = center,
                            dims = dims,
                            scale = scale)

zz[1:5,1:5]
tmp[1:5,1:5]
sum(abs(zz-tmp))

mean_vec <- Matrix::colMeans(tmp)
sd_vec <- matrixStats::colSds(tmp)
irlba_tmp2 <- irlba::irlba(tmp, nv = 5, center = NULL, scale = NULL)
irlba_tmp2$d
irlba_tmp3 <- irlba::irlba(tmp, nv = 5, center = NULL, scale = sd_vec)
irlba_tmp3$d
irlba_tmp4 <- irlba::irlba(tmp, nv = 40, center = mean_vec, scale = sd_vec)
irlba_tmp4$d
irlba_tmp5 <- irlba::irlba(tmp, nv = 5, center = mean_vec, scale = sd_vec)
irlba_tmp5$d

tmp2 <- sweep(tmp, MARGIN = 2, STATS = mean_vec, FUN = "-")
tmp2 <- sweep(tmp2, MARGIN = 2, STATS = sd_vec, FUN = "/")
rspectra_tmp <- RSpectra::svds(tmp2, k = 40)

#######



