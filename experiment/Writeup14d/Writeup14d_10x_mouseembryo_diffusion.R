rm(list=ls())
load("../../../../out/Writeup14d/Writeup14d_10x_mouseembryo_dcca.RData")

source("diffusion_distance.R")

cell_idx <- which(as.character(mbrain2@meta.data$wsnn_res.2) %in% c("17", "16", "0", "4", "2", "1", "5", "7", "12", "20"))
cell_idx <- cell_idx[-which(!mbrain2@meta.data$label_Savercat[cell_idx] %in% c("Radial glia", "Neuroblast", "Cortical or hippocampal glutamatergic"))]

set.seed(10)
res <- form_snn_graph(dcca_res, cell_idx = cell_idx)
P <- form_transition(res$snn)
diffusion_res <- extract_eigen(P, dims = 1:50, check = T)
n <- nrow(P)
dist_mat <- diffusion_distance(diffusion_res$eigenvalues, 
                               diffusion_res$right_vector)
rownames(dist_mat) <- rownames(P)
colnames(dist_mat) <- colnames(P)

############

col_palette <- scales::hue_pal()(length(unique(mbrain2@meta.data$celltype)))
col_vec <- sapply(res$cell_idx, function(x){
  col_palette[which(sort(unique(mbrain2@meta.data$celltype)) == mbrain2@meta.data$celltype[x])]
})

png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph(mbrain2[["both"]]@cell.embeddings,
           cell_idx = res$cell_idx,
           adj_mat = res$adj_mat,
           col_vec = col_vec,
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()

