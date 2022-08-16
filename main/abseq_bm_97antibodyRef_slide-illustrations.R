library(Seurat)
library(dplyr)

load("../../../out/main/abseq_bm97Ref_tcca.RData")
source("bm_97antibodyRef_colorPalette.R")

celltype_vec <- bm$ct
col_vec <- sapply(celltype_vec, function(celltype){
  col_palette[which(names(col_palette) == celltype)]
})

vec1 <- multiSVD_obj$cca_obj$score_1[,1]
vec2 <- multiSVD_obj$cca_obj$score_2[,1]

png("../../../out/figures/main/abseq_bm97Ref_cca_leading-scores.png",
    height = 1200, width = 600, units = "px", res = 300)
par(mar = c(0.5, 4, 0.5, 0.5))
plot(NA, ylim = range(c(vec1, vec2)), xlim = c(-.5, 1.5), 
     ylab = "Leading canonical score", xlab = "", bty = "n", xaxt = "n")
points(y = vec1, x = rep(0, length(vec1)) + runif(length(vec1), min = -.3, max = .3),
       pch = 16, col = col_vec, cex = 0.25)
points(y = vec2, x = rep(1, length(vec2)) + runif(length(vec1), min = -.3, max = .3), 
       pch = 16, col = col_vec, cex = 0.25)
graphics.off()

##############################

Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@scale.data[Seurat::VariableFeatures(object = bm),])
Seurat::DefaultAssay(bm) <- "AB"
mat_2 <- Matrix::t(bm[["AB"]]@scale.data)

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

set.seed(10)
multiSVD_obj2 <- tiltedCCA:::create_multiSVD(mat_1 = mat_1b, mat_2 = mat_2b,
                                             dims_1 = 1:30, dims_2 = 1:30,
                                             center_1 = T, center_2 = T,
                                             normalize_row = T,
                                             normalize_singular_value = T,
                                             recenter_1 = F, recenter_2 = F,
                                             rescale_1 = F, rescale_2 = F,
                                             scale_1 = T, scale_2 = T,
                                             verbose = 1)
multiSVD_obj2 <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj2,
                                            large_clustering_1 = as.factor(bm$RNA_snn_res.0.1), 
                                            large_clustering_2 = as.factor(bm$AB_snn_res.0.1), 
                                            num_metacells = 5000,
                                            verbose = 1)
multiSVD_obj2 <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj2,
                                          latent_k = 20,
                                          num_neigh = 15,
                                          bool_cosine = T,
                                          bool_intersect = F,
                                          min_deg = 15,
                                          verbose = 2)
multiSVD_obj3 <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj2, 
                                       fix_tilt_perc = 0.1,
                                       verbose = 1)
multiSVD_obj3 <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj3,
                                                     verbose = 1)

multiSVD_obj4 <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj2, 
                                       fix_tilt_perc = 0.9,
                                       verbose = 1)
multiSVD_obj4 <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj4,
                                                     verbose = 1)


set.seed(10)
bm[["common_tccaUp"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj3,
                                                      what = "common",
                                                      aligned_umap_assay = "rna.umap",
                                                      seurat_obj = bm,
                                                      seurat_assay = "RNA",
                                                      verbose = 1)

set.seed(10)
bm[["common_tccaDown"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj4,
                                                        what = "common",
                                                        aligned_umap_assay = "rna.umap",
                                                        seurat_obj = bm,
                                                        seurat_assay = "RNA",
                                                        verbose = 1)

bm[["common_tccaUp"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = bm[["rna.umap"]]@cell.embeddings,
  target_mat = bm[["common_tccaUp"]]@cell.embeddings
)
plot1 <- Seurat::DimPlot(bm, reduction = "common_tccaUp",
                         group.by = "ct", 
                         cols = col_palette, pt.size = 1.5)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Common UMAP 2") + ggplot2::xlab("Common UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_tcca-illustration_tiltup.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")



bm[["common_tccaDown"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = bm[["rna.umap"]]@cell.embeddings,
  target_mat = bm[["common_tccaDown"]]@cell.embeddings
)
plot1 <- Seurat::DimPlot(bm, reduction = "common_tccaDown",
                         group.by = "ct", 
                         cols = col_palette, pt.size = 1.5)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Common UMAP 2") + ggplot2::xlab("Common UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_tcca-illustration_tiltdown.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")



