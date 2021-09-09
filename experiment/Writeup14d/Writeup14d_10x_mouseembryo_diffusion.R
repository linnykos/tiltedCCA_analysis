rm(list=ls())
load("../../../../out/Writeup14d/Writeup14d_10x_mouseembryo_dcca.RData")

source("diffusion_distance.R")

# cell_idx <- cell_idx[-which(!mbrain2@meta.data$label_Savercat[cell_idx] %in% c("Radial glia", "Neuroblast", "Cortical or hippocampal glutamatergic"))]

clustering <- as.character(mbrain2@meta.data$wsnn_res.2)
selected_clusters <- c("17", "16", "0", "4", "2", "1", "5", "7", "12", "20")

set.seed(10)
res <- form_snn_graph(dcca_res, k = 5, 
                      clustering = clustering,
                      selected_clusters = selected_clusters)
P <- form_transition(res$snn,
                     lazy_param = 0.85,
                     teleport_param = 0.99)
diffusion_res <- extract_eigen(P, check = T)
n <- nrow(P)
dist_mat <- diffusion_distance(diffusion_res$eigenvalues, 
                               diffusion_res$right_vector,
                               time_vec = seq(1,5,by=1))
rownames(dist_mat) <- rownames(P)
colnames(dist_mat) <- colnames(P)

############

col_palette <- scales::hue_pal()(length(unique(clustering)))
col_vec <- sapply(rownames(res$adj_mat), function(x){
  col_palette[which(sort(unique(clustering)) == x)]
})

png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph(mbrain2[["both"]]@cell.embeddings,
           clustering = clustering,
           adj_mat = res$adj_mat,
           col_vec = col_vec,
           cex_missing = 1,
           cex_normal = 2,
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()

radial_start <- which(rownames(res$adj_mat) == "17")
cortical_end <- which.max(dist_mat[radial_start,])

# compute janky pseudotime
forward_rank <- rank(dist_mat[radial_start,])
backward_rank <- nrow(res$adj_mat) - rank(dist_mat[cortical_end,])
avg_rank <- apply(cbind(forward_rank, backward_rank), 1, mean)
avg_rank <- rank(avg_rank)

png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both_diffusion_radial.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph(mbrain2[["both"]]@cell.embeddings,
           clustering = clustering,
           adj_mat = res$adj_mat,
           feature_vec = dist_mat[radial_start,],
           zlim = range(dist_mat),
           cex_missing = 1,
           cex_normal = 2,
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()


png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both_diffusion_cortical.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph(mbrain2[["both"]]@cell.embeddings,
           clustering = clustering,
           adj_mat = res$adj_mat,
           feature_vec = dist_mat[cortical_end,],
           zlim = range(dist_mat),
           cex_missing = 1,
           cex_normal = 2,
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()


png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both_pseudotime.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph(mbrain2[["both"]]@cell.embeddings,
           clustering = clustering,
           adj_mat = res$adj_mat,
           feature_vec = avg_rank,
           cex_missing = 1,
           cex_normal = 2,
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()

