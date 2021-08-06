rm(list=ls())
load("../../../../data/Pairedtag_mousebrain_RNA-Histone/seurat_H3K4me1.RData")

Seurat::DefaultAssay(pairedtag) <- "SCT"
zz1 <- pairedtag[["SCT"]]@data[Seurat::VariableFeatures(pairedtag),]
zz1 <- Matrix::t(zz1)
mean_vec <- sparseMatrixStats::colMeans2(zz1)
tmp1a <- irlba::irlba(zz1, nv = 50)
zz1a <- multiomicCCA:::.mult_mat_vec(tmp1a$u, tmp1a$d/max(tmp1a$d))
tmp1b <- irlba::irlba(zz1, nv = 50, center = mean_vec) 
zz1b <- multiomicCCA:::.mult_mat_vec(tmp1b$u, tmp1b$d/max(tmp1b$d))
celltype <- as.factor(pairedtag@meta.data$celltype)

png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K4me1_rna_heatmap.png",
    height = 1200, width = 1800, res = 300, units = "px")
par(mar = c(4,4,4,0.5))
multiomicCCA:::plot_scores_heatmap.list(list(zz1a, zz1b), main_vec = c("RNA Uncentered", "RNA Centered"),
                                        membership_vec = celltype, log_scale = T, scaling_power = 4)
graphics.off()

Seurat::DefaultAssay(pairedtag) <- "DNA"
zz2 <- pairedtag[["DNA"]]@data[Seurat::VariableFeatures(pairedtag),]
zz2 <- Matrix::t(zz2)
mean_vec <- sparseMatrixStats::colMeans2(zz2)
tmp2a <- irlba::irlba(zz2, nv = 50)
zz2a <- multiomicCCA:::.mult_mat_vec(tmp2a$u, tmp2a$d/max(tmp2a$d))
tmp2b <- irlba::irlba(zz2, nv = 50, center = mean_vec) 
zz2b <- multiomicCCA:::.mult_mat_vec(tmp2b$u, tmp2b$d/max(tmp2b$d))

png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K4me1_histone_heatmap.png",
    height = 1200, width = 1800, res = 300, units = "px")
par(mar = c(4,4,4,0.5))
multiomicCCA:::plot_scores_heatmap.list(list(zz2a, zz2b), main_vec = c("Histone Uncentered", "Histone Centered"),
                                        membership_vec = celltype, log_scale = T, scaling_power = 5)
graphics.off()
