rm(list=ls())
load("../../../../out/Writeup14f/Writeup14f_citeseq_bm_dcca.RData")

library(Seurat); library(Signac)
library(multiomicCCA)

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)

# now investigate celltypes
set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res2, 
                                         data_1 = T, 
                                         data_2 = F,
                                         normalization_type = "cosine_itself")
set.seed(10)
adt_frnn <- multiomicCCA::construct_frnn(dcca_res2, 
                                          data_1 = F, 
                                          data_2 = T,
                                          normalization_type = "cosine_itself")
save(rna_frnn, adt_frnn, 
     file = "../../../../out/Writeup14f/Writeup14f_citeseq_bm_clisi_tmp.RData")

set.seed(10)
common_g <- multiomicCCA::combine_frnn(dcca_res2, 
                                       g_1 = rna_frnn$c_g,
                                       g_2 = adt_frnn$c_g,
                                       nn = 30,
                                       verbose = 1)
save(rna_frnn, adt_frnn, common_g,
     file = "../../../../out/Writeup14f/Writeup14f_citeseq_bm_clisi_tmp.RData")

set.seed(10)
rna_clisi <- multiomicCCA::clisi_information(common_g, rna_frnn$d_g,
                               membership_vec = factor(bm$celltype.l2),
                               max_subsample_clisi = 4000)
rna_clisi$common_clisi$clisi_mat
rna_clisi$distinct_clisi$clisi_mat

save(rna_clisi, 
     file = "../../../../out/Writeup14f/Writeup14f_citeseq_bm_clisi_tmp2.RData")

set.seed(10)
adt_clisi <- multiomicCCA::clisi_information(common_g, adt_frnn$d_g,
                                membership_vec = factor(bm$celltype.l2),
                                max_subsample_clisi = 4000)
adt_clisi$common_clisi$clisi_mat
adt_clisi$distinct_clisi$clisi_mat

save(rna_clisi, adt_clisi,
     file = "../../../../out/Writeup14f/Writeup14f_citeseq_bm_clisi_tmp2.RData")

########################
load("../../../../out/Writeup14f/Writeup14f_citeseq_bm_dcca.RData")
load("../../../../out/Writeup14f/Writeup14f_citeseq_bm_clisi_tmp2.RData")

mat_c <- rna_clisi$common_clisi$clisi_mat
mat_d1 <- rna_clisi$distinct_clisi$clisi_mat
mat_d2 <- adt_clisi$distinct_clisi$clisi_mat
color_vec <- scales::hue_pal()(nrow(mat_c))
max(c(mat_c, mat_d1, mat_d2))
max_val <- 1
k <- nrow(mat_c)
idx_list <- vector("list", 4)
idx_list[[1]] <- 1:k
idx_list[[2]] <- which(sort(unique(bm$celltype.l2)) %in% c("Prog_B 1", "Prog_B 2", "Prog_DC", "Prog_Mk", "Prog_RBC", "HSC", "LMPP", "GMP"))
idx_list[[3]] <- which(sort(unique(bm$celltype.l2)) %in% c("CD4 Memory", "CD4 Naive", "CD8 Effector_1", "CD8 Effector_2", "CD8 Memory_1", "CD8 Memory_2", "CD8 Naive"))
idx_list[[4]] <- which(sort(unique(bm$celltype.l2)) %in% "CD4 Naive")
names(idx_list) <- c("All", "Progenitors", "Developed_cells", "cd4")
  
for(kk in 1:length(idx_list)){
  png(paste0("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_clisi_", names(idx_list)[kk], ".png"), 
      height = 1300, width = 2000, res = 300, units = "px")
  par(mfrow = c(1,2), mar = c(2.5, 0.5, 2, 3))
  plot(NA, xlim = c(-max_val, 0), ylim = c(0, max_val),
       xlab = "",
       ylab = "", bty = "n", yaxt = "n", xaxt = "n", asp = T)
  axis(side = 1, at = seq(-1, 0, by = 0.2), labels = format(seq(1,0,by=-0.2), nsmall = 1),
       line = -1.5)
  axis(side = 4, at = seq(0, 1, by = 0.2), col = 2, col.axis = 2)
  graphics::mtext("Common enrichment", side = 4, line = 2.5, col = 2)
  graphics::mtext("Distinct 1 enrichment", side = 1, line = 1)
  lines(c(-100, 100), c(100, -100), col = 2, lty = 2, lwd = 2)
  
  for(i in idx_list[[kk]]){
    points(x = -mat_d1[i,i], y = mat_c[i,i], pch = 16, cex = 4.5, col = color_vec[i])
    points(x = -mat_d1[i,i], y = mat_c[i,i], pch = 16, cex = 3.5, col = "white")

    c_vec <- mat_c[,i]; c_vec <- c_vec/sum(c_vec)
    c_df <- data.frame(value = c_vec, color = color_vec, idx = 1:length(c_vec))
    c_df <- c_df[c(i, c(1:nrow(c_df))[-i]),]
    c_idx <- c(1, order(c_df$value[-1], decreasing = T)[1:2]+1)
    c_vec_final <- c(c_df$value[c_idx], sum(c_df$value[-c_idx]))
    c_col_final <- c(c_df$color[c_idx], "gray")
    d_vec <- mat_d1[,i]; d_vec <- d_vec/sum(d_vec)
    d_df <- data.frame(value = d_vec, color = color_vec, idx = 1:length(d_vec))
    d_df <- d_df[c(i, c(1:nrow(d_df))[-i]),]
    d_idx <- c(1, order(d_df$value[-1], decreasing = T)[1:2]+1)
    d_vec_final <- c(d_df$value[d_idx], sum(d_df$value[-d_idx]))
    d_col_final <- c(d_df$color[d_idx], "gray")
    pie_custom(x = c(c_vec_final, rev(d_vec_final)),
               border = c(rep(2, length(c_vec_final)), rep(1, length(d_vec_final))),
               offset = c(-mat_d1[i,i], mat_c[i,i]),
               radius = 0.045, col = c(c_col_final, rev(d_col_final)), lwd = 1)
  }
  
  par(mar = c(2.5, 3, 2, 0.5))
  plot(NA, xlim = c(0, max_val), ylim = c(0, max_val), xlab = "",
       ylab = "", bty = "n", asp = T, yaxt = "n", xaxt = "n")
  graphics::mtext("Distinct 2 enrichment", side = 1, line = 1)
  axis(side = 1, at = seq(0, 1, by = 0.2),
       line = -1.5)
  axis(side = 2, at = seq(0, 1, by = 0.2), col = 2, col.axis = 2)
  lines(c(-100, 100), c(-100, 100), col = 2, lty = 2, lwd = 2)
  for(i in idx_list[[kk]]){
    points(x = mat_d2[i,i], y = mat_c[i,i], pch = 16, cex = 4.5, col = color_vec[i])
    points(x = mat_d2[i,i], y = mat_c[i,i], pch = 16, cex = 3.5, col = "white")

    c_vec <- mat_c[,i]; c_vec <- c_vec/sum(c_vec)
    c_df <- data.frame(value = c_vec, color = color_vec, idx = 1:length(c_vec))
    c_df <- c_df[c(i, c(1:nrow(c_df))[-i]),]
    c_idx <- c(1, order(c_df$value[-1], decreasing = T)[1:2]+1)
    c_vec_final <- c(c_df$value[c_idx], sum(c_df$value[-c_idx]))
    c_col_final <- c(c_df$color[c_idx], "gray")
    d_vec <- mat_d2[,i]; d_vec <- d_vec/sum(d_vec)
    d_df <- data.frame(value = d_vec, color = color_vec, idx = 1:length(d_vec))
    d_df <- d_df[c(i, c(1:nrow(d_df))[-i]),]
    d_idx <- c(1, order(d_df$value[-1], decreasing = T)[1:2]+1)
    d_vec_final <- c(d_df$value[d_idx], sum(d_df$value[-d_idx]))
    d_col_final <- c(d_df$color[d_idx], "gray")
    pie_custom(x = c(d_vec_final, rev(c_vec_final)),
               offset = c(mat_d2[i,i], mat_c[i,i]),
               border = c(rep(1, length(d_vec_final)), rep(2, length(c_vec_final))),
               radius = 0.045, col = c(d_col_final, rev(c_col_final)), lwd = 1)
  }
  graphics.off()
}

#######

kk <- 4
png(paste0("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_clisi_", names(idx_list)[kk], "2.png"), 
    height = 1300, width = 2000, res = 300, units = "px")
par(mfrow = c(1,2), mar = c(2.5, 0.5, 2, 3))
plot(NA, xlim = c(-max_val, 0), ylim = c(0, max_val),
     xlab = "",
     ylab = "", bty = "n", yaxt = "n", xaxt = "n", asp = T)
axis(side = 1, at = seq(-1, 0, by = 0.2), labels = format(seq(1,0,by=-0.2), nsmall = 1),
     line = -1.5)
axis(side = 4, at = seq(0, 1, by = 0.2), col = 2, col.axis = 2)
graphics::mtext("Common enrichment", side = 4, line = 2.5, col = 2)
graphics::mtext("Distinct 1 enrichment", side = 1, line = 1)
lines(c(-100, 100), c(100, -100), col = 2, lty = 2, lwd = 2)

for(i in idx_list[[kk]]){
  points(x = -mat_d1[i,i], y = mat_c[i,i], pch = 16, cex = 4.5, col = color_vec[i])
  
  c_vec <- mat_c[,i]; c_vec <- c_vec/sum(c_vec)
  c_df <- data.frame(value = c_vec, color = color_vec, idx = 1:length(c_vec))
  c_df <- c_df[c(i, c(1:nrow(c_df))[-i]),]
  c_idx <- c(1, order(c_df$value[-1], decreasing = T)[1:2]+1)
  c_vec_final <- c(c_df$value[c_idx], sum(c_df$value[-c_idx]))
  c_col_final <- c(c_df$color[c_idx], "gray")
  d_vec <- mat_d1[,i]; d_vec <- d_vec/sum(d_vec)
  d_df <- data.frame(value = d_vec, color = color_vec, idx = 1:length(d_vec))
  d_df <- d_df[c(i, c(1:nrow(d_df))[-i]),]
  d_idx <- c(1, order(d_df$value[-1], decreasing = T)[1:2]+1)
  d_vec_final <- c(d_df$value[d_idx], sum(d_df$value[-d_idx]))
  d_col_final <- c(d_df$color[d_idx], "gray")
}
par(mar = c(2.5, 3, 2, 0.5))
plot(NA, xlim = c(0, max_val), ylim = c(0, max_val), xlab = "",
     ylab = "", bty = "n", asp = T, yaxt = "n", xaxt = "n")
graphics::mtext("Distinct 2 enrichment", side = 1, line = 1)
axis(side = 1, at = seq(0, 1, by = 0.2),
     line = -1.5)
axis(side = 2, at = seq(0, 1, by = 0.2), col = 2, col.axis = 2)
lines(c(-100, 100), c(-100, 100), col = 2, lty = 2, lwd = 2)
for(i in idx_list[[kk]]){
  points(x = mat_d2[i,i], y = mat_c[i,i], pch = 16, cex = 4.5, col = color_vec[i])
  
  c_vec <- mat_c[,i]; c_vec <- c_vec/sum(c_vec)
  c_df <- data.frame(value = c_vec, color = color_vec, idx = 1:length(c_vec))
  c_df <- c_df[c(i, c(1:nrow(c_df))[-i]),]
  c_idx <- c(1, order(c_df$value[-1], decreasing = T)[1:2]+1)
  c_vec_final <- c(c_df$value[c_idx], sum(c_df$value[-c_idx]))
  c_col_final <- c(c_df$color[c_idx], "gray")
  d_vec <- mat_d2[,i]; d_vec <- d_vec/sum(d_vec)
  d_df <- data.frame(value = d_vec, color = color_vec, idx = 1:length(d_vec))
  d_df <- d_df[c(i, c(1:nrow(d_df))[-i]),]
  d_idx <- c(1, order(d_df$value[-1], decreasing = T)[1:2]+1)
  d_vec_final <- c(d_df$value[d_idx], sum(d_df$value[-d_idx]))
  d_col_final <- c(d_df$color[d_idx], "gray")
}
graphics.off()

#########################################

png(paste0("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_clisi_", names(idx_list)[4], "_zoom.png"), 
    height = 1500, width = 1500, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(NA, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), xlab = "",
     ylab = "", bty = "n", asp = T, yaxt = "n", xaxt = "n")
i <- idx_list[[4]][1]
points(x = 0, y = 0, pch = 16, cex = 55, col = color_vec[i])
points(x = 0, y = 0, pch = 16, cex = 45, col = "white")

c_vec <- mat_c[,i]; c_vec <- c_vec/sum(c_vec)
c_df <- data.frame(value = c_vec, color = color_vec, idx = 1:length(c_vec))
c_df <- c_df[c(i, c(1:nrow(c_df))[-i]),]
c_idx <- c(1, order(c_df$value[-1], decreasing = T)[1:2]+1)
c_vec_final <- c(c_df$value[c_idx], sum(c_df$value[-c_idx]))
c_col_final <- c(c_df$color[c_idx], "gray")
d_vec <- mat_d1[,i]; d_vec <- d_vec/sum(d_vec)
d_df <- data.frame(value = d_vec, color = color_vec, idx = 1:length(d_vec))
d_df <- d_df[c(i, c(1:nrow(d_df))[-i]),]
d_idx <- c(1, order(d_df$value[-1], decreasing = T)[1:2]+1)
d_vec_final <- c(d_df$value[d_idx], sum(d_df$value[-d_idx]))
d_col_final <- c(d_df$color[d_idx], "gray")
pie_custom(x = c(c_vec_final, rev(d_vec_final)),
           border = c(rep(2, length(c_vec_final)), rep(1, length(d_vec_final))),
           offset = c(0, 0),
           radius = 1, col = c(d_col_final, rev(c_col_final)), lwd = 4.5)
graphics.off()

c_col_final2 <- c_col_final
c_col_final2[4] <- rgb(0.5, 0.5, 0.5)
png(paste0("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_clisi_", names(idx_list)[4], "_zoom2.png"), 
    height = 1500, width = 1500, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(NA, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), xlab = "",
     ylab = "", bty = "n", asp = T, yaxt = "n", xaxt = "n")
i <- idx_list[[4]][1]
points(x = 0, y = 0, pch = 16, cex = 55, col = color_vec[i])
points(x = 0, y = 0, pch = 16, cex = 45, col = "white")

c_vec <- mat_c[,i]; c_vec <- c_vec/sum(c_vec)
c_df <- data.frame(value = c_vec, color = color_vec, idx = 1:length(c_vec))
c_df <- c_df[c(i, c(1:nrow(c_df))[-i]),]
c_idx <- c(1, order(c_df$value[-1], decreasing = T)[1:2]+1)
c_vec_final <- c(c_df$value[c_idx], sum(c_df$value[-c_idx]))
c_col_final <- c(c_df$color[c_idx], "gray")
d_vec <- mat_d1[,i]; d_vec <- d_vec/sum(d_vec)
d_df <- data.frame(value = d_vec, color = color_vec, idx = 1:length(d_vec))
d_df <- d_df[c(i, c(1:nrow(d_df))[-i]),]
d_idx <- c(1, order(d_df$value[-1], decreasing = T)[1:2]+1)
d_vec_final <- c(d_df$value[d_idx], sum(d_df$value[-d_idx]))
d_col_final <- c(d_df$color[d_idx], "gray")
pie_custom(x = c(c_vec_final, rev(d_vec_final)),
           border = c(rep(2, length(c_vec_final)), rep(rgb(0.5, 0.5, 0.5, 0.5), length(d_vec_final))),
           offset = c(0, 0),
           radius = 1, col = c(d_col_final, rev(paste0(c_col_final2, "4D"))), lwd = 4.5)
graphics.off()

#########################################

plot1 <- Seurat::VlnPlot(bm, features = "RNA.weight", 
                         group.by = 'celltype.l2', 
                         cols = scales::hue_pal()(length(unique(bm$celltype.l2))),
                         sort = TRUE, pt.size = 0.1) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_wnn.png"),
                plot1, device = "png", width = 8, height = 4, units = "in")
