rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)

n <- ncol(bm)
svd_1 <- irlba::irlba(bm[["RNA"]]@scale.data, 30)
embedding_1 <- multiomicCCA:::.mult_mat_vec(svd_1$v, svd_1$d/svd_1$d[1]*sqrt(n))
svd_2 <- svd(bm[["ADT"]]@scale.data)
embedding_2 <- multiomicCCA:::.mult_mat_vec(svd_2$v, svd_2$d/svd_2$d[1]*sqrt(n))
embedding_all <- cbind(embedding_1, embedding_2)
rownames(embedding_all) <- colnames(bm)

embedding_all <- scale(embedding_all, center = T, scale = F)
pca_res <- svd(embedding_all)
consensus_mat <- multiomicCCA:::.mult_mat_vec(pca_res$u[,1:20], pca_res$d[1:20])

set.seed(10)
consensus_umap <- Seurat::RunUMAP(consensus_mat, 
                                  metric = "euclidean",
                                  reduction.key = "umapConsensusPCA_")
bm[["consensus.umap"]] <- Seurat::CreateDimReducObject(consensus_umap@cell.embeddings)

u_mat1 <- bm[["rna.umap"]]@cell.embeddings
u_mat2 <- bm[["consensus.umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(bm@meta.data)
colnames(tmp) <- c("umapConsensusPCA_1", "umapConsensusPCA_2")
bm[["consensus.umap"]]@cell.embeddings <- tmp


plot1 <- Seurat::DimPlot(bm, reduction = "consensus.umap", 
                         group.by = "celltype.l2",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nConsensus PCA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_consensuspca.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

###########################

vec1 <- embedding_1[,1]; vec1 <- scale(vec1, center = T, scale = F)
vec2 <- embedding_2[,1]; vec2 <- scale(vec2, center = T, scale = F)

color_df <- data.frame(celltype = sort(unique(bm@meta.data$celltype.l2)),
                       color = scales::hue_pal()(length(unique(bm@meta.data$celltype.l2))))
col_vec <- rep(NA, length(vec2))
uniq_celltypes <- sort(unique(bm@meta.data$celltype.l2))
for(celltype in uniq_celltypes){
  col_vec[which(bm@meta.data$celltype.l2 == celltype)] <- color_df[which(color_df$celltype == celltype),"color"]
}
col_vec_trans <- sapply(col_vec, function(x){
  paste0(x, "1A")
})

png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_example_pcs.png",
    height = 1200, width = 600, units = "px", res = 300)
par(mar = c(0.5, 4, 0.5, 0.5))
plot(NA, ylim = c(-2,2), xlim = c(-.5, 1.5), 
     ylab = "Leading PCs", xlab = "", bty = "n", xaxt = "n")
points(y = vec1, x = rep(0, length(vec1)) + runif(length(vec1), min = -.3, max = .3),
       pch = 16, col = col_vec, cex = 0.25)
points(y = vec2, x = rep(1, length(vec2)) + runif(length(vec1), min = -.3, max = .3), 
       pch = 16, col = col_vec, cex = 0.25)
graphics.off()

xlim <- c(-2, 6); ylim <- c(-5, 3)
png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_consensus_leadingpc.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(5,5,0.5,0.5))
plot(NA, asp = T, 
     xlim = xlim, ylim = ylim, bty = "n",
     xlab = expression(Leading ~ PC ~ of ~ Modality ~ 1: ~ sigma^(1)==175),
     ylab = expression(Leading ~ PC ~ of ~ Modality ~ 2: ~ sigma^(2)==175))
lines(c(-1e5,1e5), rep(0,2),lty = 2)
lines(rep(0,2), c(-1e5,1e5), lty = 2)
points(vec1, vec2, col = col_vec, pch = 16,)
graphics.off()

pca_res <- stats::prcomp(cbind(vec1, vec2))
n <- length(vec1)
set.seed(10)
jitter_vec <- runif(n, min = -.1, max = .1)
jitter_2d <- jitter_vec %*% t(pca_res$rotation[,2])
projection_vec <- cbind(vec1, vec2) %*% pca_res$rotation[,1] %*% t(pca_res$rotation[,1])
projection_vec <- projection_vec + jitter_2d
png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_consensus_leadingpc_projected.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(5,5,0.5,0.5))
plot(NA, asp = T, 
     xlim = xlim, ylim = ylim, bty = "n",
     xlab = expression(Leading ~ PC ~ of ~ Modality ~ 1: ~ sigma^(1)==175),
     ylab = expression(Leading ~ PC ~ of ~ Modality ~ 2: ~ sigma^(2)==175))
lines(c(-1e5,1e5), rep(0,2),lty = 2)
lines(rep(0,2), c(-1e5,1e5), lty = 2)
points(vec1, vec2, col = col_vec_trans, pch = 16)
points(projection_vec, pch = 16, col = "white", cex = 2)
points(projection_vec, pch = 16, col = col_vec)
graphics.off()