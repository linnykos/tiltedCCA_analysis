rm(list=ls())
load("../../../../out/Writeup14d/Writeup14d_10x_mouseembryo_dcca.RData")

source("diffusion_distance.R")

cell_idx <- which(as.character(mbrain2@meta.data$wsnn_res.2) %in% c("17", "16", "0", "4", "2", "1", "5", "7", "12", "20"))
cell_idx <- cell_idx[-which(!mbrain2@meta.data$label_Savercat[cell_idx] %in% c("Radial glia", "Neuroblast", "Cortical or hippocampal glutamatergic"))]

set.seed(10)
res <- form_snn_graph(dcca_res, k = 5, cell_idx = cell_idx)
P <- form_transition(res$snn,
                     lazy_param = 0.85,
                     teleport_param = 0.99)
diffusion_res <- extract_eigen(P, dims = 1:100, check = T)
n <- nrow(P)
dist_mat <- diffusion_distance(diffusion_res$eigenvalues, 
                               diffusion_res$right_vector,
                               time_vec = seq(1,1000,by=25))
rownames(dist_mat) <- rownames(P)
colnames(dist_mat) <- colnames(P)

############

col_palette <- scales::hue_pal()(length(unique(mbrain2@meta.data$celltype)))
col_vec <- sapply(res$cell_idx, function(x){
  col_palette[which(sort(unique(mbrain2@meta.data$celltype)) == mbrain2@meta.data$celltype[x])]
})

png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both.png",
    height = 2500, width = 2500, units = "px", res = 300)
plot_graph(mbrain2[["both"]]@cell.embeddings,
           cell_idx = res$cell_idx,
           adj_mat = res$adj_mat,
           col_vec = col_vec,
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()

radial_idx <- intersect(which(mbrain2@meta.data$celltype == "Radial glia"),
                        which(as.character(mbrain2@meta.data$wsnn_res.2) == "17"))
radial_idx2 <- which(rownames(res$adj_mat) %in% rownames(mbrain2@meta.data[radial_idx,]))
cortical_idx <- intersect(which(mbrain2@meta.data$celltype == "Cortical or hippocampal glutamatergic"),
                          which(as.character(mbrain2@meta.data$wsnn_res.2) == "20"))
cortical_idx2 <- which(rownames(res$adj_mat) %in% rownames(mbrain2@meta.data[cortical_idx,]))

starting_radial <- radial_idx2[which.max(apply(dist_mat[radial_idx2,cortical_idx2], 1, median))]
ending_cortical <- cortical_idx2[which.max(apply(dist_mat[cortical_idx2,radial_idx2], 1, median))]

# compute janky pseudotime
forward_rank <- rank(dist_mat[starting_radial,])
backward_rank <- n - rank(dist_mat[ending_cortical,])
avg_rank <- apply(cbind(forward_rank, backward_rank), 1, mean)
avg_rank <- rank(avg_rank)

png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both_diffusion_radial.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph(mbrain2[["both"]]@cell.embeddings,
           cell_idx = res$cell_idx,
           adj_mat = res$adj_mat,
           feature_vec = dist_mat[starting_radial,],
           zlim = range(dist_mat),
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()


png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both_diffusion_cortical.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph(mbrain2[["both"]]@cell.embeddings,
           cell_idx = res$cell_idx,
           adj_mat = res$adj_mat,
           feature_vec = dist_mat[ending_cortical,],
           zlim = range(dist_mat),
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()


png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both_pseudotime.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph(mbrain2[["both"]]@cell.embeddings,
           cell_idx = res$cell_idx,
           adj_mat = res$adj_mat,
           feature_vec = avg_rank,
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()

