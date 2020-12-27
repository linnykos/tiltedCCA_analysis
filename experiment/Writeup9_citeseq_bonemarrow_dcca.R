rm(list=ls())
set.seed(10)

date_of_run <- Sys.time()
session_info <- sessionInfo()

library(Seurat)
load("../../out/Writeup9_citeseq_bonemarrow_dimred_all.RData")

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)

library(multiomicCCA)

# choosing a rank
svd_res_1 <- RSpectra::svds(mat_1, k = 25)
svd_res_2 <- svd(mat_2)
-diff(svd_res_1$d)/svd_res_1$d[-1]; -diff(svd_res_2$d)/svd_res_2$d[-1]

# doing D-CCA
K <- 8
dcca_res <- multiomicCCA::dcca(mat_1, mat_2, rank_1 = K, rank_2 = K, rank_12 = K, 
                               enforce_rank = T, verbose = T)

save.image("../../out/Writeup9_citeseq_bonemarrow_dcca.RData")

##########################

# plot ranks
png("../../out/Writeup9_citeseq_bonemarrow_dcca_rank.png", height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,2))
graphics::plot(svd_res_1$d, xlab = "Rank", ylab = "Singular value", main = "SVD for RNA", pch = 16)
graphics::lines(rep(K+.5, 2), c(-1e5, 1e5), col = "red", lty = 2, lwd = 2)
graphics::plot(svd_res_2$d, xlab = "Rank", ylab = "Singular value", main = "SVD for Protein", pch = 16)
graphics::lines(rep(K+.5, 2), c(-1e5, 1e5), col = "red", lty = 2, lwd = 2)
graphics.off()

# plot ranks of decompositions
common_svd_1 <- RSpectra::svds(dcca_res$common_mat_1, k = K)$d
distinct_svd_1 <- RSpectra::svds(dcca_res$distinct_mat_1, k = K)$d
common_svd_2 <- RSpectra::svds(dcca_res$common_mat_2, k = K)$d
distinct_svd_2 <- RSpectra::svds(dcca_res$distinct_mat_2, k = K)$d

png("../../out/Writeup9_citeseq_bonemarrow_dcca_explained_variability.png", height = 900, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3))
graphics::plot(common_svd_1, ylim = range(c(0, common_svd_1, distinct_svd_1)), xlab = "Rank", ylab = "Singular value", 
               main = paste0("Explained variability for RNA\n(Total explained by common: ",
                             round(100*sum(common_svd_1)/sum(c(common_svd_1,distinct_svd_1)), 2), "%)"), pch = 16, cex = 2)
graphics::points(distinct_svd_1, pch = 16, cex = 2, col = "red")
legend("topright", c("Common", "Distinct"), bty="n", fill=c("black", "red"))
graphics::plot(common_svd_2, ylim = range(c(0, common_svd_2, distinct_svd_2)), xlab = "Rank", ylab = "Singular value", 
               main = paste0("Explained variability for Protein\n(Total explained by common: ",
               round(100*sum(common_svd_2)/sum(c(common_svd_2,distinct_svd_2)), 2), "%)"), pch = 16, cex = 2)
graphics::points(distinct_svd_2, pch = 16, cex = 2, col = "red")
legend("topright", c("Common", "Distinct"), bty="n", fill=c("black", "red"))
graphics::plot(dcca_res$cca_obj, ylim = c(0,1), xlab = "Rank", ylab = "Canonical correlation", 
               main = "CCA objective", pch = 16, cex = 2)
graphics.off()

#################

# first plot the usual UMAPs suggested by the Seurat vingette
set.seed(10)
bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

png("../../out/Writeup9_citeseq_bonemarrow_seurat_wnn_umap.png", height = 1500, width = 1500, units = "px", res = 300)
# par(mfrow = c(2,2))
plot1 <- Seurat::DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("WNN UMAP (Seurat)")
graphics.off()

png("../../out/Writeup9_citeseq_bonemarrow_seurat_rna_umap.png", height = 1500, width = 1500, units = "px", res = 300)
# par(mfrow = c(2,2))
plot1 <- Seurat::DimPlot(bm, reduction = 'rna.umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("RNA UMAP (Vanilla)")
graphics.off()

png("../../out/Writeup9_citeseq_bonemarrow_seurat_protein_umap.png", height = 1500, width = 1500, units = "px", res = 300)
# par(mfrow = c(2,2))
plot1 <- Seurat::DimPlot(bm, reduction = 'adt.umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("Protein UMAP (Vanilla)")
graphics.off()

##################


# plot embeddings
# we're going to essentially try all possible variants:
# - based on only the common factors
# - based on RNA's weighted common SVD and distinctive weighted SVD
# - based on Protein's weighted common SVD and distinctive weighted SVD
# - based on combining all 4
# first extract all the necessary 
bm2 <- bm
common_fact <- dcca_res$common_factors
rownames(common_fact) <- colnames(bm)
set.seed(10)
zz <- Seurat::RunUMAP(common_fact, reduction.key = 'DCCA_common_')
bm2[["common_factor"]] <- zz

png("../../out/Writeup9_citeseq_bonemarrow_dcca_common_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm2, reduction = 'common_factor', group.by = 'celltype.l2', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("Common UMAP (D-CCA)")
graphics.off()

svd_tmp <- RSpectra::svds(dcca_res$common_mat_1+dcca_res$distinct_mat_1, k = K)
tmp <- multiomicCCA:::.mult_mat_vec(svd_tmp$u, svd_tmp$d)
rownames(tmp) <- colnames(bm)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCARNA_')
bm2[["dcca_rna_umap"]] <- zz
png("../../out/Writeup9_citeseq_bonemarrow_dcca_rna_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm2, reduction = 'dcca_rna_umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("RNA-only UMAP (D-CCA)")
graphics.off()

svd_tmp <- RSpectra::svds(dcca_res$common_mat_2+dcca_res$distinct_mat_2, k = K)
tmp <- multiomicCCA:::.mult_mat_vec(svd_tmp$u, svd_tmp$d)
rownames(tmp) <- colnames(bm)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCAProtein_')
bm2[["dcca_protein_umap"]] <- zz
png("../../out/Writeup9_citeseq_bonemarrow_dcca_protein_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm2, reduction = 'dcca_protein_umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("Protein-only UMAP (D-CCA)")
graphics.off()

# now do the throw-everything-together UMAP
svd_list <- vector("list", length = 4)
svd_list[[1]] <- RSpectra::svds(dcca_res$common_mat_1, k = K)
svd_list[[2]] <- RSpectra::svds(dcca_res$common_mat_2, k = K)
svd_list[[3]] <- RSpectra::svds(dcca_res$distinct_mat_1, k = K)
svd_list[[4]] <- RSpectra::svds(dcca_res$distinct_mat_2, k = K)

tmp <- do.call(cbind, lapply(svd_list, function(svd_tmp){
  multiomicCCA:::.mult_mat_vec(svd_tmp$u, svd_tmp$d)
}))
rownames(tmp) <- colnames(bm)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'DCCAEverything_')
bm2[["dcca_everything_umap"]] <- zz
png("../../out/Writeup9_citeseq_bonemarrow_dcca_everything_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm2, reduction = 'dcca_everything_umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("Everything UMAP (D-CCA)")
graphics.off()

##############################################

# plot WNN or D-CCA RNA-weights
png("../../out/Writeup9_citeseq_bonemarrow_seurat_weight_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::FeaturePlot(bm, features = "RNA.weight", 
                  reduction = 'wnn.umap')
plot1 + ggplot2::ggtitle("Seurat WNN's RNA weight")
graphics.off()

png("../../out/Writeup9_citeseq_bonemarrow_seurat_weight.png", height = 1500, width = 2500, units = "px", res = 300)
Seurat::VlnPlot(bm, features = "RNA.weight", group.by = 'celltype.l2', sort = TRUE, pt.size = 0.1) +
  NoLegend() + ggplot2::ggtitle("Seurat WNN's RNA weight")
graphics.off()

# now compute our weights
.l2norm <- function(x){sqrt(sum(x^2))}
dcca_weights <- sapply(1:nrow(dcca_res$common_mat_1), function(i){
  if(i %% floor(nrow(dcca_res$common_mat_1)/10) == 0) cat('*')
  rna_weight <- .l2norm(dcca_res$distinct_mat_1[i,])/(.l2norm(dcca_res$common_mat_1[i,])+.l2norm(dcca_res$distinct_mat_1[i,]))
  protein_weight <- .l2norm(dcca_res$distinct_mat_2[i,])/(.l2norm(dcca_res$common_mat_2[i,])+.l2norm(dcca_res$distinct_mat_2[i,]))
  exp(rna_weight)/(exp(rna_weight)+exp(protein_weight))
})
bm2@meta.data[["dcca_weight"]] <- dcca_weights
png("../../out/Writeup9_citeseq_bonemarrow_dcca_weight_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::FeaturePlot(bm2, features = "dcca_weight", 
                             reduction = 'dcca_everything_umap')
plot1 + ggplot2::ggtitle("D-CCA's RNA weight")
graphics.off()

png("../../out/Writeup9_citeseq_bonemarrow_dcca_weight.png", height = 1500, width = 2500, units = "px", res = 300)
Seurat::VlnPlot(bm2, features = "dcca_weight", group.by = 'celltype.l2', sort = TRUE, pt.size = 0.1) +
  NoLegend() + ggplot2::ggtitle("D-CCA's RNA weight")
graphics.off()

##################################

svd_1 <- RSpectra::svds(mat_1, k = K)
svd_2 <- RSpectra::svds(mat_2, k = K)
mat_total <- cbind(mat_1/svd_1$d[1], mat_2/svd_2$d[1])

pca_total <- RSpectra::svds(mat_total, k = 2*K)
tmp <- multiomicCCA:::.mult_mat_vec(pca_total$u, pca_total$d)
rownames(tmp) <- colnames(bm)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'PCAFused_')
bm2[["pca_fused_umap"]] <- zz
png("../../out/Writeup9_citeseq_bonemarrow_pca_fused_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm2, reduction = 'pca_fused_umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("PCA fused UMAP")
graphics.off()


tmp <- cbind(multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d), multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d))
rownames(tmp) <- colnames(bm)
set.seed(10)
zz <- Seurat::RunUMAP(tmp, reduction.key = 'PCAConcat_')
bm2[["pca_concate_umap"]] <- zz
png("../../out/Writeup9_citeseq_bonemarrow_pca_concate_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm2, reduction = 'pca_concate_umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + NoLegend()
plot1 + ggplot2::ggtitle("PCA-concatenated UMAP")
graphics.off()
