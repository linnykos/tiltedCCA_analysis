rm(list=ls())
library(Seurat)
load("../../../out/main/citeseq_bm25_tcca.RData")
source("bm_25antibody_colorPalette.R")

# consensus pca
Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@data[Seurat::VariableFeatures(object = bm),])
Seurat::DefaultAssay(bm) <- "ADT"
mat_2 <- Matrix::t(bm[["ADT"]]@data[Seurat::VariableFeatures(object = bm),])

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
consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = mat_1b, mat_2 = mat_2b,
                                           dims_1 = 1:30, dims_2 = 1:18,
                                           dims_consensus = 1:30,
                                           center_1 = T, center_2 = T,
                                           recenter_1 = F, recenter_2 = F,
                                           rescale_1 = F, rescale_2 = F,
                                           scale_1 = T, scale_2 = T,
                                           verbose = 1)
consensus_dimred <- consensus_pca$dimred_consensus

############################

membership_vec <- factor(bm$celltype.l2)

cell_enrichment_res <- tiltedCCA:::postprocess_cell_enrichment(
  input_obj = multiSVD_obj, 
  membership_vec = membership_vec, 
  max_subsample = 1000,
  verbose = 1
)
cell_enrichment_res$enrichment_common$df

consensus_enrichment <- tiltedCCA:::postprocess_cell_enrichment(
  input_obj = consensus_dimred, 
  membership_vec = membership_vec, 
  num_neigh = multiSVD_obj$param$snn_num_neigh,
  bool_cosine = multiSVD_obj$param$snn_bool_cosine,
  bool_intersect = multiSVD_obj$param$snn_bool_intersect,
  max_subsample = 1000,
  min_deg = multiSVD_obj$param$snn_min_deg,
  verbose = 1
)
consensus_enrichment$enrichment$df

#####################

x_vec <- cell_enrichment_res$enrichment_common$df[,"value"]
names(x_vec) <- cell_enrichment_res$enrichment_common$df[,"celltype"]
y_vec <- consensus_enrichment$enrichment$df[,"value"]
names(y_vec) <- consensus_enrichment$enrichment$df[,"celltype"]

png("../../../out/figures/main/bm_25antibody_consensusPCA-tCCA_enrichment.png",
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(4, 4, 0.5, 0.5))
plot(NA, xlim = c(0,1), ylim = c(0,1),
     main = "", xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n", asp = T)
for(x in seq(0,1,by=0.1)){
  lines(rep(x,2), c(-10,10), col = "gray", lty = 3)
}
for(y in seq(0,1,by=0.1)){
  lines(c(-10,10), rep(y,2), col = "gray", lty = 3)
}
lines(c(-10,10), c(-10,10), col = 2, lty = 2, lwd = 2)

col_vec <- col_palette[names(x_vec)]
points(x_vec, y_vec, col = 1, cex = 4, pch = 16)
points(x_vec, y_vec, col = "white", cex = 3, pch = 16)
points(x_vec, y_vec, col = col_vec, cex = 2.5, pch = 16)

axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)

graphics.off()
