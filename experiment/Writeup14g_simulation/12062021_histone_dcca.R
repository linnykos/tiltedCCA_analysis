rm(list=ls())
load("../../out/simulation/Writeup14g_data.RData")

mat_1 <- seurat_obj2[["pca"]]@cell.embeddings
mat_2 <- seurat_obj2[["lsi"]]@cell.embeddings
rank_1 <- 20; rank_2 <- 20; nn <- 30
set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, 
                                      dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      num_neigh = nn, 
                                      metacell_clustering_1 = NA,
                                      metacell_clustering_2 = NA,
                                      fix_tilt_perc = T, verbose = T)
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, rank_c = 20)

set.seed(10)
dcca_common_umap <- Seurat::RunUMAP(cbind(dcca_decomp$common_mat_1, dcca_decomp$common_mat_2), 
                                    metric = "euclidean",
                                    reduction.key = "dccaCommon_")
seurat_obj2[["dcca_common"]] <- Seurat::CreateDimReducObject(dcca_common_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap1"]]@cell.embeddings
u_mat2 <- seurat_obj2[["dcca_common"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("dccaCommon_1", "dccaCommon_2")
seurat_obj2[["dcca_common"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "dcca_common", 
                         cols = color_vec, 
                         group.by = "celltype_custom",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: D-CCA (Common)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/dcca_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

color_vec2 <- rep("gray", length(color_vec))
names(color_vec2) <- names(color_vec)
plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "dcca_common", 
                         cols = color_vec2)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: D-CCA (Common)"))
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/dcca_common_none.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

##################################


vec1 <- dcca_res$score_1[,1]
vec2 <- dcca_res$score_2[,1]
col_vec <- rep(NA, length(vec1))
uniq_celltypes <- sort(unique(seurat_obj2@meta.data$celltype))
for(celltype in uniq_celltypes){
  col_vec[which(seurat_obj2@meta.data$celltype == celltype)] <- color_df[which(color_df$celltype == celltype),"color"]
}

png("../../out/simulation/Writeup14e_simulation/example_scores.png",
    height = 1200, width = 600, units = "px", res = 300)
par(mar = c(0.5, 4, 0.5, 0.5))
plot(NA, ylim = range(c(vec1, vec2)), xlim = c(-.5, 1.5), 
     ylab = "Leading canonical score", xlab = "", bty = "n", xaxt = "n")
points(y = vec1, x = rep(0, length(vec1)) + runif(length(vec1), min = -.3, max = .3),
       pch = 16, col = col_vec, cex = 0.25)
points(y = vec2, x = rep(1, length(vec2)) + runif(length(vec1), min = -.3, max = .3), 
       pch = 16, col = col_vec, cex = 0.25)
graphics.off()
