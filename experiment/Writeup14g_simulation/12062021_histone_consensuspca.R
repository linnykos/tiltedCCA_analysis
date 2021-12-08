rm(list=ls())
load("../../out/simulation/Writeup14g_data.RData")

mat1 <- seurat_obj2[["pca"]]@cell.embeddings
mat2 <- seurat_obj2[["lsi"]]@cell.embeddings
vec1 <- seurat_obj2[["pca"]]@cell.embeddings[,1]
vec2 <- seurat_obj2[["lsi"]]@cell.embeddings[,1]

pca_res <- stats::prcomp(cbind(vec1, vec2))

xlim <- c(-0.8, 1) # quantile(c(vec1, vec2), probs = c(0.05, 0.95))

col_vec <- rep(NA, length(vec2))
uniq_celltypes <- sort(unique(seurat_obj2@meta.data$celltype))
for(celltype in uniq_celltypes){
  col_vec[which(seurat_obj2@meta.data$celltype == celltype)] <- color_df[which(color_df$celltype == celltype),"color"]
}
col_vec_trans <- sapply(col_vec, function(x){
  paste0(x, "1A")
})

png("../../out/figures/Writeup14g_simulation/consensus_leadingpc.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(5,5,0.5,0.5))
plot(NA, asp = T, 
     xlim = xlim, ylim = xlim, bty = "n",
     xlab = expression(Leading ~ PC ~ of ~ Modality ~ 1: ~ sigma^(1)==20),
     ylab = expression(Leading ~ PC ~ of ~ Modality ~ 2: ~ sigma^(2)==20))
lines(c(-1e5,1e5), rep(0,2),lty = 2)
lines(rep(0,2), c(-1e5,1e5), lty = 2)
points(vec1, vec2, col = col_vec, pch = 16,)
graphics.off()


n <- length(vec1)
set.seed(10)
jitter_vec <- runif(n, min = -.1, max = .1)
jitter_2d <- jitter_vec %*% t(pca_res$rotation[,2])
projection_vec <- cbind(vec1, vec2) %*% pca_res$rotation[,1] %*% t(pca_res$rotation[,1])
projection_vec <- projection_vec + jitter_2d
png("../../out/figures/Writeup14g_simulation/consensus_leadingpc_projected.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(5,5,0.5,0.5))
plot(NA, asp = T, 
     xlim = xlim, ylim = xlim, bty = "n",
     xlab = expression(Leading ~ PC ~ of ~ Modality ~ 1: ~ sigma^(1)==20),
     ylab = expression(Leading ~ PC ~ of ~ Modality ~ 2: ~ sigma^(2)==20))
lines(c(-1e5,1e5), rep(0,2),lty = 2)
lines(rep(0,2), c(-1e5,1e5), lty = 2)
points(vec1, vec2, col = col_vec_trans, pch = 16)
points(projection_vec, pch = 16, col = "white", cex = 2)
points(projection_vec, pch = 16, col = col_vec)
graphics.off()


png("../../out/figures/Writeup14g_simulation/example_pcs.png",
    height = 1200, width = 600, units = "px", res = 300)
par(mar = c(0.5, 4, 0.5, 0.5))
plot(NA, ylim = c(-0.8,1), xlim = c(-.5, 1.5), 
     ylab = "Leading PCs", xlab = "", bty = "n", xaxt = "n")
points(y = vec1, x = rep(0, length(vec1)) + runif(length(vec1), min = -.3, max = .3),
       pch = 16, col = col_vec, cex = 0.25)
points(y = vec2, x = rep(1, length(vec2)) + runif(length(vec1), min = -.3, max = .3), 
       pch = 16, col = col_vec, cex = 0.25)
graphics.off()

#################################

n <- ncol(seurat_obj2)
svd_1 <- svd(seurat_obj2[["pca"]]@cell.embeddings)
embedding_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
svd_2 <- svd(seurat_obj2[["lsi"]]@cell.embeddings)
embedding_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))
embedding_all <- cbind(embedding_1, embedding_2)
rownames(embedding_all) <- colnames(seurat_obj2)
pca_res <- stats::prcomp(embedding_all, center = TRUE, scale. = FALSE)
consensus_mat <- pca_res$x[,1:40] # Hm... I think this should be 20?

set.seed(10)
consensus_umap <- Seurat::RunUMAP(consensus_mat, 
                                  metric = "euclidean",
                                  reduction.key = "umapConsensusPCA_")
seurat_obj2[["consensus.umap"]] <- Seurat::CreateDimReducObject(consensus_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap1"]]@cell.embeddings
u_mat2 <- seurat_obj2[["consensus.umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("umapConsensusPCA_1", "umapConsensusPCA_2")
seurat_obj2[["consensus.umap"]]@cell.embeddings <- tmp


plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "consensus.umap", 
                         cols = color_vec, 
                         group.by = "celltype_custom",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Consensus PCA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/consensuspca.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


