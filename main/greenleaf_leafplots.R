rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_greenleaf_tcca_RNA-geneActivity2.RData")
multiSVD_obj2 <- multiSVD_obj
greenleaf2 <- greenleaf
load("../../../out/main/10x_greenleaf_tcca_RNA-geneActivity.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

main_genes <- c("ASCL1", "NEUROD1", "NTRK2")
secondary_genes <- c("EOMES", "NHLH1", "PAX6", "ZNF521",
                     "HES4", "WNT11", "MEIS2", "KCNH8",
                     "CNTN1", "GRM1", "NEFL", "SEMA3E")
var_features <- c(main_genes, secondary_genes)
all(var_features %in% colnames(multiSVD_obj$common_mat_1))
all(var_features %in% colnames(multiSVD_obj2$common_mat_1))

# extract relevant matrix for Seurat Gene activities
cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_common <- multiSVD_obj$common_mat_1[cell_idx[order_vec],var_features]
atac_common <- multiSVD_obj$common_mat_2[cell_idx[order_vec],paste0("ATAC-", var_features)]

######

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_mat <- multiSVD_obj$common_mat_1 + multiSVD_obj$distinct_mat_1
rna_mat <- rna_mat[cell_idx[order_vec],var_features] 

mat <- Matrix::t(greenleaf[["customGAct"]]@data[rownames(multiSVD_obj$svd_2$v),])
mat <- scale(mat)
set.seed(10)
svd_atac <- irlba::irlba(mat, nv = 50)
atac_mat <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_atac$u[,2:50], svd_atac$d[2:50]), svd_atac$v[,2:50])
rownames(atac_mat) <- rownames(mat)
colnames(atac_mat) <- colnames(mat)
atac_mat <- atac_mat[cell_idx[order_vec],paste0("ATAC-",var_features)]
all(colnames(atac_mat) == paste0("ATAC-", colnames(rna_mat)))

##############################

# extract relevant matrix for CICERO gene activities
cell_idx <- intersect(which(greenleaf2$Lineage1 == 1), which(!is.na(greenleaf2$pseudotime)))
pseudotime_vec <- greenleaf2$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_common2 <- multiSVD_obj2$common_mat_1[cell_idx[order_vec],var_features]
atac_common2 <- multiSVD_obj2$common_mat_2[cell_idx[order_vec],paste0("ATAC-", var_features)]

########

cell_idx <- intersect(which(greenleaf2$Lineage1 == 1), which(!is.na(greenleaf2$pseudotime)))
pseudotime_vec <- greenleaf2$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_mat2 <- multiSVD_obj2$common_mat_1 + multiSVD_obj2$distinct_mat_1
rna_mat2 <- rna_mat2[cell_idx[order_vec],var_features] 

mat <- Matrix::t(greenleaf2[["geneActivity"]]@data[rownames(multiSVD_obj2$svd_2$v),])
mat <- scale(mat)
set.seed(10)
svd_atac <- irlba::irlba(mat, nv = 50)
atac_mat2 <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_atac$u[,2:50], svd_atac$d[2:50]), svd_atac$v[,2:50])
rownames(atac_mat2) <- rownames(mat)
colnames(atac_mat2) <- colnames(mat)
atac_mat2 <- atac_mat2[cell_idx[order_vec],paste0("ATAC-",var_features)]
all(colnames(atac_mat2) == paste0("ATAC-", colnames(rna_mat2)))

#################################

construction_smoothed_matrix <- function(mat){
  n <- nrow(mat)
  p <- ncol(mat)
  
  pred_mat <- sapply(1:p, function(j){
    print(j)
    
    tmp_df <- data.frame(y = mat[,j], x = 1:n)
    reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
    x_vec <- 1:n
    y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
    y_vec$Estimation[,"Pred"]
  })
  
  colnames(pred_mat) <- colnames(mat)
  pred_mat
}

pred_rna_common <- construction_smoothed_matrix(rna_common)
pred_atac_common <- construction_smoothed_matrix(atac_common)
pred_rna_mat <- construction_smoothed_matrix(rna_mat)
pred_atac_mat <- construction_smoothed_matrix(atac_mat)

pred_rna_common2 <- construction_smoothed_matrix(rna_common2)
pred_atac_common2 <- construction_smoothed_matrix(atac_common2)
pred_rna_mat2 <- construction_smoothed_matrix(rna_mat2)
pred_atac_mat2 <- construction_smoothed_matrix(atac_mat2)

README <- "rna_common is using Seurat gene activity, while rna_common2 is using CICERO. And common is using tilted-cca, and mat is original"

save(rna_common, atac_common, rna_mat, atac_mat,
     rna_common2, atac_common2, rna_mat2, atac_mat2,
     pred_rna_common, pred_atac_common, pred_rna_mat, pred_atac_mat,
     pred_rna_common2, pred_atac_common2, pred_rna_mat2, pred_atac_mat2,
     date_of_run, session_info, README,
     file = "../../../out/main/10x_greenleaf_leafplot_data.RData")

#############################################
#############################################
#############################################
#############################################

rm(list=ls())
load("../../out/main/10x_greenleaf_leafplot_data.RData")

n <- nrow(rna_mat)
color_vec <- rev(grDevices::colorRampPalette(c(rgb(191, 74, 223, maxColorValue = 255),
                                               rgb(0.8, 0.8, 0.8),
                                               rgb(239, 158, 88, maxColorValue = 255)))(n))
example_color <- rgb(1,1,1,0.5)
color_vec_trans <- sapply(color_vec, function(x){
  paste0(x, substr(example_color, start = 8, stop = 10))
})

# plot(1:n, col= color_vec_trans, pch = 16)

main_genes <- c("ASCL1", "NEUROD1", "NTRK2")
secondary_genes <- c("EOMES", "NHLH1", "PAX6", "ZNF521",
                     "HES4", "WNT11", "MEIS2", "KCNH8",
                     "CNTN1", "GRM1", "NEFL", "SEMA3E")

for(gene in main_genes){
  png(paste0("../../out/figures/main/10x_greenleaf_leafplot_mainplot_seuratGActivity_common_", gene, ".png"),
      height = 2000, width = 2000, units = "px", res = 500)
  par(mar = c(2.5,2.5,0.5,0.5))
  set.seed(10)
  idx <- sample(1:n)
  
  x_vec1 <- atac_common[,paste0("ATAC-", gene)]
  x_vec2 <- pred_atac_common[,paste0("ATAC-", gene)]
  y_vec1 <- rna_common[,gene]
  y_vec2 <- pred_rna_common[,gene]
  
  x_mean <- mean(x_vec2); x_vec1 <- x_vec1- x_mean; x_vec2 <- x_vec2 - x_mean
  x_sd <- sd(x_vec2); x_vec1 <- x_vec1/x_sd; x_vec2 <- x_vec2/x_sd
  y_mean <- mean(y_vec2); y_vec1 <- y_vec1- y_mean; y_vec2 <- y_vec2 - y_mean
  y_sd <- sd(y_vec2); y_vec1 <- y_vec1/y_sd; y_vec2 <- y_vec2/y_sd
  
  plot(x = x_vec1[idx], y = y_vec1[idx],
       pch = 16, col = color_vec_trans[idx], cex = 0.75,
       xlim = quantile(x_vec1, probs = c(0.005, 0.995)),
       ylim = quantile(y_vec1, probs = c(0.005, 0.995)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n")
  axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
  axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
  
  points(x = x_vec2[idx], 
         y = y_vec2[idx], 
         pch = 16, col = "white", cex = 2)
  points(x = rev(x_vec2), 
         y = rev(y_vec2), 
         pch = 16, col = rev(color_vec), cex = 1)
  
  points(x = x_vec2[1], 
         y = y_vec2[1],
         pch = 16, col = "white", cex = 4)
  points(x = x_vec2[1], 
         y = y_vec2[1], 
         pch = 16, col = color_vec[1], cex = 3)
  
  x_right <- mean(x_vec1[x_vec1 >= quantile(x_vec1, probs = 0.9)])
  x_left <- mean(x_vec1[x_vec1 <= quantile(x_vec1, probs = 0.1)])
  y_top <- mean(y_vec1[y_vec1 >= quantile(y_vec1, probs = 0.9)])
  y_bot <- mean(y_vec1[y_vec1 <= quantile(y_vec1, probs = 0.1)])
  
  slope <- (y_top-y_bot)/(x_right-x_left)
  intercept <- y_top - slope*x_right
  x_left2 <- -1e4; x_right2 <- 1e4
  y_bot2 <- slope*x_left2+intercept; y_top2 <- slope*x_right2+intercept
  lines(c(x_right2,x_left2), c(y_top2, y_bot2), lwd = 2, lty = 2)
  graphics.off()
}


for(gene in main_genes){
  png(paste0("../../out/figures/main/10x_greenleaf_leafplot_mainplot_seuratGActivity_original_", gene, ".png"),
      height = 1000, width = 1000, units = "px", res = 500)
  par(mar = c(0.5,0.5,0.5,0.5))
  set.seed(10)
  idx <- sample(1:n)
  
  x_vec1 <- atac_mat[,paste0("ATAC-", gene)]
  x_vec2 <- pred_atac_mat[,paste0("ATAC-", gene)]
  y_vec1 <- rna_mat[,gene]
  y_vec2 <- pred_rna_mat[,gene]
  
  x_mean <- mean(x_vec2); x_vec1 <- x_vec1- x_mean; x_vec2 <- x_vec2 - x_mean
  x_sd <- sd(x_vec2); x_vec1 <- x_vec1/x_sd; x_vec2 <- x_vec2/x_sd
  y_mean <- mean(y_vec2); y_vec1 <- y_vec1- y_mean; y_vec2 <- y_vec2 - y_mean
  y_sd <- sd(y_vec2); y_vec1 <- y_vec1/y_sd; y_vec2 <- y_vec2/y_sd
  
  plot(x = x_vec1[idx], y = y_vec1[idx],
       pch = 16, col = color_vec_trans[idx], cex = 1,
       xlim = quantile(x_vec1, probs = c(0.005, 0.995)),
       ylim = quantile(y_vec1, probs = c(0.005, 0.995)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n")
  axis(1, labels = F, cex.axis = 2, cex.lab = 2, lwd = 2)
  axis(2, labels = F, cex.axis = 2, cex.lab = 2, lwd = 2)
  
  points(x = x_vec2[idx], 
         y = y_vec2[idx], 
         pch = 16, col = "white", cex = 3)
  points(x = rev(x_vec2), 
         y = rev(y_vec2), 
         pch = 16, col = rev(color_vec), cex = 2)
  
  points(x = x_vec2[1], 
         y = y_vec2[1],
         pch = 16, col = "white", cex = 5)
  points(x = x_vec2[1], 
         y = y_vec2[1], 
         pch = 16, col = color_vec[1], cex = 4)
  
  x_right <- mean(x_vec1[x_vec1 >= quantile(x_vec1, probs = 0.9)])
  x_left <- mean(x_vec1[x_vec1 <= quantile(x_vec1, probs = 0.1)])
  y_top <- mean(y_vec1[y_vec1 >= quantile(y_vec1, probs = 0.9)])
  y_bot <- mean(y_vec1[y_vec1 <= quantile(y_vec1, probs = 0.1)])
  
  slope <- (y_top-y_bot)/(x_right-x_left)
  intercept <- y_top - slope*x_right
  x_left2 <- -1e4; x_right2 <- 1e4
  y_bot2 <- slope*x_left2+intercept; y_top2 <- slope*x_right2+intercept
  lines(c(x_right2,x_left2), c(y_top2, y_bot2), lwd = 4, lty = 2)
  graphics.off()
}

##########

for(gene in secondary_genes){
  png(paste0("../../out/figures/main/10x_greenleaf_leafplot_mainplot_seuratGActivity_common_", gene, ".png"),
      height = 1500, width = 1500, units = "px", res = 500)
  par(mar = c(0.5,0.5,0.5,0.5))
  set.seed(10)
  idx <- sample(1:n)
  
  x_vec1 <- atac_common[,paste0("ATAC-", gene)]
  x_vec2 <- pred_atac_common[,paste0("ATAC-", gene)]
  y_vec1 <- rna_common[,gene]
  y_vec2 <- pred_rna_common[,gene]
  
  x_mean <- mean(x_vec2); x_vec1 <- x_vec1- x_mean; x_vec2 <- x_vec2 - x_mean
  x_sd <- sd(x_vec2); x_vec1 <- x_vec1/x_sd; x_vec2 <- x_vec2/x_sd
  y_mean <- mean(y_vec2); y_vec1 <- y_vec1- y_mean; y_vec2 <- y_vec2 - y_mean
  y_sd <- sd(y_vec2); y_vec1 <- y_vec1/y_sd; y_vec2 <- y_vec2/y_sd
  
  plot(x = x_vec1[idx], y = y_vec1[idx],
       pch = 16, col = color_vec_trans[idx], cex = 1.5,
       xlim = quantile(x_vec1, probs = c(0.005, 0.995)),
       ylim = quantile(y_vec1, probs = c(0.005, 0.995)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n", bty = "n")
  axis(1, labels = F, cex.axis = 2.5, cex.lab = 2.5, lwd = 2.5)
  axis(2, labels = F, cex.axis = 2.5, cex.lab = 2.5, lwd = 2.5)
  
  points(x = x_vec2[idx], 
         y = y_vec2[idx], 
         pch = 16, col = "white", cex = 4)
  points(x = rev(x_vec2), 
         y = rev(y_vec2), 
         pch = 16, col = rev(color_vec), cex = 3)
  
  points(x = x_vec2[1], 
         y = y_vec2[1],
         pch = 16, col = "white", cex = 6)
  points(x = x_vec2[1], 
         y = y_vec2[1], 
         pch = 16, col = color_vec[1], cex = 5)
  
  x_right <- mean(x_vec1[x_vec1 >= quantile(x_vec1, probs = 0.9)])
  x_left <- mean(x_vec1[x_vec1 <= quantile(x_vec1, probs = 0.1)])
  y_top <- mean(y_vec1[y_vec1 >= quantile(y_vec1, probs = 0.9)])
  y_bot <- mean(y_vec1[y_vec1 <= quantile(y_vec1, probs = 0.1)])
  
  slope <- (y_top-y_bot)/(x_right-x_left)
  intercept <- y_top - slope*x_right
  x_left2 <- -1e4; x_right2 <- 1e4
  y_bot2 <- slope*x_left2+intercept; y_top2 <- slope*x_right2+intercept
  lines(c(x_right2,x_left2), c(y_top2, y_bot2), lwd = 5, lty = 2)
  graphics.off()
}

############################################

# library(quantreg)
# library(quantreg.nonpar)
library(quantdr)

basis_bsp <- fda::create.bspline.basis(breaks=seq(0, n, length.out = 11))

rna_top_mat <- sapply(main_genes, function(gene){
  print(gene)
  
  vec <- rna_common[,gene]
  n <- length(vec)
  
  reg_res <- quantdr::llqr(x = 1:n, y = vec, tau = 0.9)
  reg_res$ll_est
})
rna_bot_mat <- sapply(main_genes, function(gene){
  print(gene)
  
  vec <- rna_common[,gene]
  n <- length(vec)
  
  reg_res <- quantdr::llqr(x = 1:n, y = vec, tau = 0.1)
  reg_res$ll_est
})
atac_top_mat <- sapply(main_genes, function(gene){
  print(gene)
  
  vec <- atac_common[,paste0("ATAC-", gene)]
  n <- length(vec)
  
  reg_res <- quantdr::llqr(x = 1:n, y = vec, tau = 0.9)
  reg_res$ll_est
})
atac_bot_mat <- sapply(main_genes, function(gene){
  print(gene)
  
  vec <- atac_common[,paste0("ATAC-", gene)]
  n <- length(vec)
  
  reg_res <- quantdr::llqr(x = 1:n, y = vec, tau = 0.1)
  reg_res$ll_est
})

for(gene in main_genes){
  png(paste0("../../out/figures/main/10x_greenleaf_timeseries_seuratGActivity_common_", gene, ".png"),
      height = 1500, width = 3500, units = "px", res = 500)
  par(mar = c(2.5,2.5,0.5,0.5))
  
  x_vec <- pred_atac_common[,paste0("ATAC-", gene)]
  y_vec <- pred_rna_common[,gene]
  x_top <- atac_top_mat[,gene]
  x_bot <- atac_bot_mat[,gene]
  y_top <- rna_top_mat[,gene]
  y_bot <- rna_bot_mat[,gene]
  
  x_mean <- mean(x_vec); x_vec <- x_vec - x_mean
  x_top <- x_top - x_mean; x_bot <- x_bot - x_mean
  x_sd <- sd(x_vec); x_vec <- x_vec/x_sd
  x_top <- x_top/x_sd; x_bot <- x_bot/x_sd
  
  y_mean <- mean(y_vec); y_vec <- y_vec - y_mean
  y_top <- y_top - y_mean; y_bot <- y_bot - y_mean
  y_sd <- sd(y_vec); y_vec <- y_vec/y_sd
  y_top <- y_top/y_sd; y_bot <- y_bot/y_sd
  
  xlim <- c(1, n)
  ylim <- range(c(x_top, x_bot, y_top, y_bot))
  
  plot(NA, xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n", bty = "n",
       xlab = "", ylab = "", main = "")
  polygon(x = c(1:n, n:1),
          y = c(y_top, rev(y_bot)),
          col = rgb(55, 100, 234, 0.4*255, maxColorValue = 255),
          border = NA)
  polygon(x = c(1:n, n:1),
          y = c(x_top, rev(x_bot)),
          col = rgb(191, 74, 223, 0.4*255, maxColorValue = 255),
          border = NA)
  
  spacing_idx <- round(seq(1, n, length.out = 500))
  lines(x = 1:n, y = y_vec,
        col = "white", lwd = 15)
  lines(x = c(1:n)[spacing_idx], y = y_vec[spacing_idx],
        col = rgb(55, 100, 234, maxColorValue = 255),
        lwd = 10)
  lines(x = 1:n, y = x_vec,
        col = "white", lwd = 15)
  lines(x = c(1:n)[spacing_idx], y = x_vec[spacing_idx],
        col = rgb(191, 74, 223, maxColorValue = 255),
        lty = 2, lwd = 10)
  
  axis(1, cex.axis = 2, cex.lab = 1.5, lwd = 3)
  axis(2, cex.axis = 2, cex.lab = 1.5, lwd = 3)
  
  graphics.off()
}

