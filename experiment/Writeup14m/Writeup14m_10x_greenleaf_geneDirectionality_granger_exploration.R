rm(list=ls())
load("../../out/experiment/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_grangerTest.RData")

n <- nrow(rna_common2)
concave_fit <- .shape_constrained_fit(rna_common2[,1],
                                      grid_vec = seq(0,n,length.out=5)[-c(1,5)],
                                      metric = 2,
                                      setting = "concave")
convex_fit <- .shape_constrained_fit(rna_common2[,1],
                                     grid_vec = seq(0,n,length.out=5)[-c(1,5)],
                                     metric = 2,
                                     setting = "convex")
increasing_fit <- .shape_constrained_fit(rna_common2[,1],
                                         grid_vec = seq(0,n,length.out=5)[-c(1,5)],
                                         metric = 2,
                                         setting = "increasing")
decreasing_fit <- .shape_constrained_fit(rna_common2[,1],
                                         grid_vec = seq(0,n,length.out=5)[-c(1,5)],
                                         metric = 2,
                                         setting = "decreasing")
vexcave_fit <- .shape_constrained_fit(rna_common2[,1],
                                      grid_vec = seq(0,n,length.out=5)[-c(1,5)],
                                      metric = 2,
                                      setting = "convex-concave")
cavevex_fit <- .shape_constrained_fit(rna_common2[,1],
                                      grid_vec = seq(0,n,length.out=5)[-c(1,5)],
                                      metric = 2,
                                      setting = "concave-convex")

res_list <- list(concave_fit, convex_fit,
                 increasing_fit, decreasing_fit,
                 vexcave_fit, cavevex_fit)
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(rna_common2[,1], pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
for(i in 1:length(res_list)){
  lines(res_list[[i]], lwd = 5, col = i)
}

#####################
#####################
#####################

j <- gene_idx_upward[1]
n <- nrow(rna_unimodal_fit_mat)
col_vec <- grDevices::colorRampPalette(c('lightgrey', 'blue'))(n)

par(mar = c(0.5,0.5,0.5,0.5))
plot(rna_unimodal_fit_mat[,j])
points(atac_unimodal_fit_mat[,j], col = 2, pch = 16)

.substantial_unimodal(threshold_perc = 0.1,
                      vec = atac_unimodal_fit_mat[,j])

plot(rna_unimodal_fit_mat[,j], atac_unimodal_fit_mat[,j],
     col = col_vec, pch = 16, asp = T)

rna_diff <- diff(rna_unimodal_fit_mat[,j])
atac_diff <- diff(atac_unimodal_fit_mat[,j])

plot(rna_diff, atac_diff)