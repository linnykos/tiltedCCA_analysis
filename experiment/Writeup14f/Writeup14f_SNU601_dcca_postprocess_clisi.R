rm(list=ls())
load("../../../../out/Writeup14f/Writeup14f_SNU601_dcca.RData")
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)

# now investigate celltypes
set.seed(10)
atac_frnn <- multiomicCCA::construct_frnn(dcca_res2, 
                                          data_1 = T, 
                                          data_2 = F,
                                          normalization_type = "signac_itself")
set.seed(10)
cna_frnn <- multiomicCCA::construct_frnn(dcca_res2, 
                                         data_1 = F, 
                                         data_2 = T,
                                         normalization_type = "cosine_itself")
set.seed(10)
common_g <- multiomicCCA::combine_frnn(dcca_res2, 
                                       g_1 = atac_frnn$c_g,
                                       g_2 = cna_frnn$c_g,
                                       nn = 30)

set.seed(10)
atac_clisi <- multiomicCCA::clisi_information(common_g, atac_frnn$d_g,
                                membership_vec = factor(SNU$clone),
                                max_subsample_clisi = 4000)
atac_clisi$common_clisi$clisi_mat
atac_clisi$distinct_clisi$clisi_mat

set.seed(10)
cna_clisi <- clisi_information(common_g, cna_frnn$d_g,
                                membership_vec = factor(SNU$clone),
                                max_subsample_clisi = 4000)
cna_clisi$common_clisi$clisi_mat
cna_clisi$distinct_clisi$clisi_mat


save(atac_clisi, cna_clisi,
     file = "../../../../out/Writeup14f/Writeup14f_SNU601_clisi_tmp.RData")
metadata <- SNU@meta.data
save(metadata,
     file = "../../../../out/Writeup14f/Writeup14f_SNU601_metadata.RData")

# ###########################
# 
# tmp_g <- multiomicCCA:::.symmetrize_sparse(atac_frnn$c_g, set_ones = T)
# idx <- which(SNU$clone == "1")
# zz <- tmp_g[idx,]
# vec <- apply(zz, 1, function(x){
#   idx2 <- which(x != 0)
#   length(intersect(idx2, idx))/length(idx2)
# })
# quantile(vec)
# table(SNU$clone)/length(SNU$clone)
# 
# i <- idx[which.max(vec)]
# atac_clisi$common_clisi$clisi_cell_mat[,which(colnames(atac_clisi$common_clisi$clisi_cell_mat)== rownames(SNU@meta.data)[i])]
# 

########################################

# for local mac
load("../../out/experiment/Writeup14f/Writeup14f_SNU601_clisi_tmp.RData")
load("../../out/experiment/Writeup14f/Writeup14f_SNU601_metadata.RData")
color_vec <- scales::hue_pal()(6)

mat_c <- atac_clisi$common_clisi$clisi_mat
mat_d1 <- atac_clisi$distinct_clisi$clisi_mat
mat_d2 <- cna_clisi$distinct_clisi$clisi_mat
k <- nrow(mat_c)

png("../../out/figures/Writeup14f/Writeup14f_SNU_clisi.png", height = 1300, 
    width = 2000, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(4, 0.5, 2, 3))
max_val <- max(c(mat_c, mat_d1, mat_d2))
plot(NA, xlim = c(-max_val, 0), ylim = c(0, max_val), 
     xlab = "Distinct 1 enrichment",
     ylab = "", bty = "n", yaxt = "n", xaxt = "n", asp = T)
axis(side = 1, at = seq(-0.5, 0, by = 0.1), labels = as.character(seq(0.5, 0, by = -0.1)))
axis(side = 4, at = seq(0, 0.5, by = 0.1), col = 2, col.axis = 2)
graphics::mtext("Common enrichment", side = 4, line = 2.5, col = 2)
for(i in 1:k){
  points(x = -mat_d1[i,i], y = mat_c[i,i], pch = 16, cex = 4.5, col = color_vec[i])
  points(x = -mat_d1[i,i], y = mat_c[i,i], pch = 16, cex = 3.5, col = "white")
  
  c_vec <- mat_c[,i]; c_vec <- c_vec/sum(c_vec)
  d_vec <- mat_d1[,i]; d_vec <- d_vec/sum(d_vec)
  pie_custom(x = c(c_vec, rev(d_vec)),
             offset = c(-mat_d1[i,i], mat_c[i,i]),
             radius = 0.025, col = c(color_vec, rev(color_vec)), lwd = 1)
}
lines(c(-100, 100), c(100, -100), col = 2, lty = 2, lwd = 2)

par(mar = c(4, 3, 2, 0.5))
plot(NA, xlim = c(0, max_val), ylim = c(0, max_val), xlab = "Distinct 2 enrichment",
     ylab = "", bty = "n", asp = T, yaxt = "n")
axis(side = 2, at = seq(0, 0.5, by = 0.1), col = 2, col.axis = 2)
for(i in 1:k){
  points(x = mat_d2[i,i], y = mat_c[i,i], pch = 16, cex = 4.5, col = color_vec[i])
  points(x = mat_d2[i,i], y = mat_c[i,i], pch = 16, cex = 3.5, col = "white")
  
  c_vec <- mat_c[,i]; c_vec <- c_vec/sum(c_vec)
  d_vec <- mat_d2[,i]; d_vec <- d_vec/sum(d_vec)
  pie_custom(x = c(d_vec, rev(c_vec)),
             offset = c(mat_d2[i,i], mat_c[i,i]),
             border = c(rep(1, length(d_vec)), rep(2, length(c_vec))),
             radius = 0.025, col = c(color_vec, rev(color_vec)), lwd = 1)
}
lines(c(-100, 100), c(-100, 100), col = 2, lty = 2, lwd = 2)

graphics.off()

###############################

idx <- which(metadata$clone == "5")
row_idx <- 5
mat_c <- atac_clisi$common_clisi$clisi_cell_mat[row_idx,idx]
mat_d1 <- atac_clisi$distinct_clisi$clisi_cell_mat[row_idx,idx]
mat_d2 <- cna_clisi$distinct_clisi$clisi_cell_mat[row_idx,idx]

z_quantiles <- c(0, 0.5, 0.65, 0.8, 0.9, 0.95, 1)
z_levels <- quantile(f1$z, probs = z_quantiles)
density_colors <- colorRampPalette(c("white", "gray"))(6)
f1 <- MASS::kde2d(x = -mat_d1, y = mat_c, n = 100)
# min_val <- min(f1$z); max_val <- max(f1$z); tol <- 1e-2
# z_levels <- pretty(c(min_val, max_val), n = 11)
image(f1, col = density_colors, xlim = c(-1,0), ylim = c(0,1),
      breaks = z_levels)
graphics::contour(f1, add = T, 
                  drawlabels = F, levels = z_levels)
points(x = -mat_d1, y = mat_c, col = rgb(0.5,0.5,0.5,0.5))

par(mfrow = c(1,2), mar = c(4, 0.5, 2, 3))
max_val <- max(c(mat_c, mat_d1, mat_d2))
plot(NA, xlim = c(-max_val, 0), ylim = c(0, max_val), 
     xlab = "",
     ylab = "", bty = "n", yaxt = "n", xaxt = "n", asp = T)
axis(side = 1, at = seq(-1, 0, by = 0.2), labels = format(seq(1,0,by=-0.2), nsmall = 1),
     line = -1.5)
axis(side = 4, at = seq(0, 1, by = 0.2), col = 2, col.axis = 2)
graphics::mtext("Common enrichment", side = 4, line = 2.5, col = 2)
graphics::mtext("Distinct 1 enrichment", side = 1, line = 1)
lines(c(-100, 100), c(100, -100), col = 2, lty = 2, lwd = 2)
for(i in 1:length(mat_c)){
  points(x = -mat_d1[i], y = mat_c[i], pch = 16, cex = 1, col = rgb(0.5,0.5,0.5,0.5))
}

par(mar = c(4, 3, 2, 0.5))
plot(NA, xlim = c(0, max_val), ylim = c(0, max_val), xlab = "",
     ylab = "", bty = "n", asp = T, yaxt = "n", xaxt = "n")
graphics::mtext("Distinct 2 enrichment", side = 1, line = 1)
axis(side = 1, at = seq(0, 1, by = 0.2),
     line = -1.5)
axis(side = 2, at = seq(0, 1, by = 0.2), col = 2, col.axis = 2)
lines(c(-100, 100), c(-100, 100), col = 2, lty = 2, lwd = 2)
for(i in 1:length(mat_c)){
  points(x = mat_d2[i], y = mat_c[i], pch = 16, cex = 1, col = rgb(0.5,0.5,0.5,0.5))
}

