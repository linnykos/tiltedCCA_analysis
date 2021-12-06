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
cna_clisi <- multiomicCCA::clisi_information(common_g, cna_frnn$d_g,
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

# see https://www.r-graph-gallery.com/100-high-density-scatterplot-with-binning.html
idx <- which(metadata$clone == "5")
row_idx <- 5
mat_c <- atac_clisi$common_clisi$clisi_cell_mat[row_idx,idx]
mat_d1 <- atac_clisi$distinct_clisi$clisi_cell_mat[row_idx,idx]
mat_d2 <- cna_clisi$distinct_clisi$clisi_cell_mat[row_idx,idx]

png("../../out/figures/Writeup14f/Writeup14f_SNU_clisi_clone5_part1.png", height = 1300, 
    width = 1500, res = 300, units = "px")
par(mar = c(4,4,0.5,0.5))
bin <- hexbin::hexbin(-mat_d1, mat_c)
my_colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,'Spectral')))
plot(bin, main="" , colramp = my_colors, trans = function(x){sqrt(x)},
     inv = function(x){x^2},
     xlab = "Distinct 1 enrichment",
     ylab = "Common enrichment")
graphics.off()
png("../../out/figures/Writeup14f/Writeup14f_SNU_clisi_clone5_part2.png", height = 1300, 
    width = 1500, res = 300, units = "px")
par(mar = c(4,4,0.5,0.5))
bin <- hexbin::hexbin(mat_d2, mat_c)
my_colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,'Spectral')))
plot(bin, main="" , colramp = my_colors, trans = function(x){sqrt(x)},
     inv = function(x){x^2},
     xlab = "Distinct 2 enrichment",
     ylab = "Common enrichment")
graphics.off()

#########################

idx <- which(metadata$clone == "6")
row_idx <- 6
mat_c <- atac_clisi$common_clisi$clisi_cell_mat[row_idx,idx]
mat_d1 <- atac_clisi$distinct_clisi$clisi_cell_mat[row_idx,idx]
mat_d2 <- cna_clisi$distinct_clisi$clisi_cell_mat[row_idx,idx]

png("../../out/figures/Writeup14f/Writeup14f_SNU_clisi_clone6_part1.png", height = 1300, 
    width = 1500, res = 300, units = "px")
par(mar = c(4,4,0.5,0.5))
bin <- hexbin::hexbin(-mat_d1, mat_c)
my_colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,'Spectral')))
plot(bin, main="" , colramp = my_colors, trans = function(x){sqrt(x)},
     inv = function(x){x^2},
     xlab = "Distinct 1 enrichment",
     ylab = "Common enrichment")
graphics.off()
png("../../out/figures/Writeup14f/Writeup14f_SNU_clisi_clone6_part2.png", height = 1300, 
    width = 1500, res = 300, units = "px")
par(mar = c(4,4,0.5,0.5))
bin <- hexbin::hexbin(mat_d2, mat_c)
my_colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,'Spectral')))
plot(bin, main="" , colramp = my_colors, trans = function(x){sqrt(x)},
     inv = function(x){x^2},
     xlab = "Distinct 2 enrichment",
     ylab = "Common enrichment")
graphics.off()

##########################

load("../../out/experiment/Writeup14f/Writeup14f_SNU601_umap.RData")
zz <- SNU2[["dcca_distinct1"]]@cell.embeddings
cell_idx <- which(metadata$clone == "5")
xrange <- c(0, 2.5)
yrange <- c(-4, -1.5)
plot(zz[,1], zz[,2], asp = T, col = rgb(0.5,0.5,0.5,0.5), pch = 16)
points(zz[cell_idx,1], zz[cell_idx,2], col = 1, pch = 16)
cell_idx2 <- which(sapply(cell_idx, function(i){
  bool1 <- zz[i,1] >= xrange[1]
  bool2 <- zz[i,1] <= xrange[2]
  bool3 <- zz[i,2] >= yrange[1]
  bool4 <- zz[i,2] <= yrange[2]
  
  bool1 & bool2 & bool3 & bool4
}))
points(zz[cell_idx[cell_idx2],1], zz[cell_idx[cell_idx2],2], col = 2, pch = 16)
cell_names <- rownames(SNU2@meta.data)[cell_idx[cell_idx2]]

mat_c <- atac_clisi$common_clisi$clisi_cell_mat[5,cell_names]
mat_d1 <- atac_clisi$distinct_clisi$clisi_cell_mat[5,cell_names]
mat_d2 <- cna_clisi$distinct_clisi$clisi_cell_mat[5,cell_names]

png("../../out/figures/Writeup14f/Writeup14f_SNU_clisi_clone5_part1_sub.png", height = 1300, 
    width = 1500, res = 300, units = "px")
par(mar = c(4,4,0.5,0.5))
bin <- hexbin::hexbin(-mat_d1, mat_c)
my_colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,'Spectral')))
plot(bin, main="" , colramp = my_colors, trans = function(x){sqrt(x)},
     inv = function(x){x^2},
     xlab = "Distinct 1 enrichment",
     ylab = "Common enrichment")
graphics.off()


zz <- SNU2[["dcca_distinct2"]]@cell.embeddings
cell_idx <- which(metadata$clone == "6")
xrange <- c(-5,-3)
yrange <- c(-5, -3)
plot(zz[,1], zz[,2], asp = T, col = rgb(0.5,0.5,0.5,0.5), pch = 16)
points(zz[cell_idx,1], zz[cell_idx,2], col = 1, pch = 16)
cell_idx2 <- which(sapply(cell_idx, function(i){
  bool1 <- zz[i,1] >= xrange[1]
  bool2 <- zz[i,1] <= xrange[2]
  bool3 <- zz[i,2] >= yrange[1]
  bool4 <- zz[i,2] <= yrange[2]
  
  bool1 & bool2 & bool3 & bool4
}))
points(zz[cell_idx[cell_idx2],1], zz[cell_idx[cell_idx2],2], col = 2, pch = 16)
cell_names <- rownames(SNU2@meta.data)[cell_idx[cell_idx2]]

mat_c <- atac_clisi$common_clisi$clisi_cell_mat[6,cell_names]
mat_d1 <- atac_clisi$distinct_clisi$clisi_cell_mat[6,cell_names]
mat_d2 <- cna_clisi$distinct_clisi$clisi_cell_mat[6,cell_names]

png("../../out/figures/Writeup14f/Writeup14f_SNU_clisi_clone6_part2_sub.png", height = 1300, 
    width = 1500, res = 300, units = "px")
par(mar = c(4,4,0.5,0.5))
bin <- hexbin::hexbin(mat_d2, mat_c)
my_colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,'Spectral')))
plot(bin, main="" , colramp = my_colors, trans = function(x){sqrt(x)},
     inv = function(x){x^2},
     xlab = "Distinct 2 enrichment",
     ylab = "Common enrichment")
graphics.off()
