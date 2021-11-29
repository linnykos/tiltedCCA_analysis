rm(list=ls())
load("../../../../out/Writeup14f/Writeup14f_SNU601_dcca.RData")
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)

library(mclust); library(Seurat); library(Signac)

summary_mat <- compute_variable_summary(mat = dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2,
                                        common_mat = dcca_decomp$common_mat_2,
                                        metacell_clustering = factor(SNU$clone),
                                        verbose = 1)

col_vec <- rep(1, ncol(dcca_decomp$common_mat_2))
col_vec[grep("theta", colnames(dcca_decomp$common_mat_2))] <- 2

png("../../../../out/figures/Writeup14f/Writeup14f_SNU601_cna_exploration.png", 
    height = 1500, width = 1500, res = 300, units = "px")
plot(summary_mat[,2], summary_mat[,1], pch = 16,
     col = col_vec, xlab = "Separability (KL divergence)",
     ylab = "Alignment w/ common space (R^2)",
     main = "SNU601 (CNA):\n(Theta: red, Rho: black)")
graphics.off()

grep("-3:",  colnames(dcca_decomp$common_mat_2))
grep("-20:",  colnames(dcca_decomp$common_mat_2))
idx_3_rho <- c(7,8,9)
idx_3_theta <- c(71,72,73)
idx_20_rho <- c(59,60)
idx_20_theta <- c(123,124)

png("../../../../out/figures/Writeup14f/Writeup14f_SNU601_cna_exploration_specific.png", 
    height = 1500, width = 1500, res = 300, units = "px")
plot(summary_mat[,2], summary_mat[,1], pch = 16,
     col = "gray", xlab = "Separability (KL divergence)",
     ylab = "Alignment w/ common space (R^2)",
     main = "Chr3: red, Chr20: green\nRho: circle, Theta: triangle")
points(summary_mat[c(idx_3_rho, idx_3_theta, idx_20_rho, idx_20_theta),2], 
       summary_mat[c(idx_3_rho, idx_3_theta, idx_20_rho, idx_20_theta),1],
       pch = c(rep(16,3),rep(17,2), rep(16,3),rep(17,2)),
       col = c(rep(2, 5), rep(3, 5)),
       cex = 2)
graphics.off()

rownames(summary_mat)[order(summary_mat[,2], decreasing = T)[1:10]]

#####################
load("../../../../out/Writeup14f/Writeup14f_SNU_preprocessed.RData")
Seurat::DefaultAssay(SNU) <- "cna"
mat_2 <- Matrix::t(SNU[["cna"]]@data)

j_vec <- c(8,60)
uniq_clone <- sort(unique(SNU$clone))
color_vec <- scales::hue_pal()(length(uniq_clone))

for(j in 1:2){
    png(paste0("../../../../out/figures/Writeup14f/Writeup14f_SNU601_cna_exploration_specific_",
        rownames(summary_mat)[j_vec[j]], ".png"), 
        height = 2000, width = 1500, 
        res = 300, units = "px")
    par(mfrow = c(3,2), mar = c(4,4,4,0.5))
    xlim <- c(min(mat_2[,j_vec[j]]), quantile(mat_2[,j_vec[j]], prob = 0.95))
    for(i in 1:length(uniq_clone)){
        cell_idx <- intersect(
            intersect(which(SNU$clone == uniq_clone[i]),
                      which(mat_2[,j_vec[j]] >= xlim[1])),
            which(mat_2[,j_vec[j]] <= xlim[2]))
        hist(mat_2[cell_idx,j_vec[j]], col = color_vec[i], xlim = xlim,
             main = paste0("Clone ", uniq_clone[i]), breaks = 25)
    }
    graphics.off()
}
