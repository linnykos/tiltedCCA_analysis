rm(list=ls())
library(Seurat); library(Signac)

load("../../../../out/Writeup14p/Writeup14p_10x_greenleaf_tcca.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)

##################################

mentioned_genes <- c("ALDH2", "APOE", "AQP4", "ASCL1", "BHLHE22",
                     "BHLHE40", "C16orf89", "CAV2", "DLX2", "DOK5",
                     "DUSP1", "EOMES", "ETV4", "FOS", "FOXJ1", "GLI3",
                     "HAS2", "HES1", "HES4", "HSPA1A", "HSPA1B",
                     "ID3", "IGFBP7", "JUN", "KIF1A", "LIMCH1",
                     "MBP", "MEF2C", "NEUROD1", "NEUROD2", "NEUROD4",
                     "NEUROD6", "NEUROG1", "NEUROG2", "NFIA", "NFIB", "NFIC",
                     "NHLH1", "NR2F1", "PAX6", "RFX4", "RUNX1",
                     "OLIG1", "OLIG2", "SOX2", "SOX3", "SOX6",
                     "SOX9", "SOX10", "SOX21", "SPARCL1", "SNCB", "TBX",
                     "TNC", "TOP2A", "TRB1", "WNT11")
mentioned_genes <- mentioned_genes[mentioned_genes %in% colnames(rna_common)]

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_common <- rna_common[cell_idx[order_vec],mentioned_genes]
atac_pred <- atac_pred[cell_idx[order_vec],mentioned_genes]

p <- ncol(rna_common); n <- nrow(rna_common)
rna_common_smoothed <- sapply(1:p, function(j){
  print(j)
  tmp_df <- data.frame(y = rna_common[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
atac_pred_smoothed <- sapply(1:p, function(j){
  print(j)
  tmp_df <- data.frame(y = atac_pred[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})

colnames(rna_common_smoothed) <- colnames(rna_common)
colnames(atac_pred_smoothed) <- colnames(atac_pred)

##########

n <- nrow(rna_common)
for(k in 1:2){
  genes <- mentioned_genes[((k-1)*30+1):min((k*30), length(mentioned_genes))]
  
  png(paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_greenleaf_timeseries_npregfast_mentioned_",
             k, ".png"),
      height = 4500, width = 3000, units = "px", res = 300)
  par(mfrow = c(6,5), mar = c(4,4,4,0.5))
  for(gene in genes){
    set.seed(10)
    x_vec <- c(1:n, 1:n)
    y_vec <- c(rna_common[,gene], atac_pred[,gene])
    col <- c(rep(rgb(0.5, 0.5, 0.5, 0.2), n), rep(rgb(0.8, 0, 0, 0.2), n))
    idx <- sample(1:length(x_vec))
    plot(x = x_vec[idx], y = y_vec[idx], 
         pch = 16, col = col[idx], main = gene)
    lines(rna_common_smoothed[,gene], col = "white", lwd = 4)
    lines(rna_common_smoothed[,gene], col = 1, lwd = 2)
    
    lines(atac_pred_smoothed[,gene], col = "white", lwd = 4)
    lines(atac_pred_smoothed[,gene], col = 2, lwd = 2)
  }
  graphics.off()
}

################################################
################################################
################################################

rna_all <- multiSVD_obj$common_mat_1 +  multiSVD_obj$distinct_mat_1
rna_all <- rna_all[cell_idx[order_vec],mentioned_genes]
rna_all_smoothed <- sapply(1:p, function(j){
  print(j)
  tmp_df <- data.frame(y = rna_all[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(rna_all_smoothed) <- colnames(rna_all)

n <- nrow(rna_common)
for(k in 1:2){
  genes <- mentioned_genes[((k-1)*30+1):min((k*30), length(mentioned_genes))]
  
  png(paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_greenleaf_timeseries_npregfast_mentioned_",
             k, "_v2.png"),
      height = 4500, width = 3000, units = "px", res = 300)
  par(mfrow = c(6,5), mar = c(4,4,4,0.5))
  for(gene in genes){
    set.seed(10)
    x_vec <- c(1:n, 1:n)
    y_vec <- c(rna_all[,gene], atac_pred[,gene])
    col <- c(rep(rgb(0.5, 0.5, 0.5, 0.2), n), rep(rgb(0.8, 0, 0, 0.2), n))
    idx <- sample(1:length(x_vec))
    plot(x = x_vec[idx], y = y_vec[idx], 
         pch = 16, col = col[idx], main = gene)
    lines(rna_all_smoothed[,gene], col = "white", lwd = 4)
    lines(rna_all_smoothed[,gene], col = 1, lwd = 2)
    
    lines(atac_pred_smoothed[,gene], col = "white", lwd = 4)
    lines(atac_pred_smoothed[,gene], col = 2, lwd = 2)
  }
  graphics.off()
}

#############################

svd_1 <- tiltedCCA:::.svd_safe(mat = multiSVD_obj$common_mat_1,
                               check_stability = T,
                               K = 49, # positive integer
                               mean_vec = NULL, # boolean, NULL or vector
                               rescale = F, # boolean
                               scale_max = NULL, # NULL or positive integer
                               sd_vec = NULL)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred2 <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
atac_pred2 <- atac_pred2[cell_idx[order_vec],mentioned_genes]

atac_pred2_smoothed <- sapply(1:p, function(j){
  print(j)
  tmp_df <- data.frame(y = atac_pred2[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(atac_pred2_smoothed) <- colnames(atac_pred2)

n <- nrow(rna_common)
for(k in 1:2){
  genes <- mentioned_genes[((k-1)*30+1):min((k*30), length(mentioned_genes))]
  
  png(paste0("../../../../out/figures/Writeup14p/Writeup14p_10x_greenleaf_timeseries_npregfast_mentioned_",
             k, "_v3.png"),
      height = 4500, width = 3000, units = "px", res = 300)
  par(mfrow = c(6,5), mar = c(4,4,4,0.5))
  for(gene in genes){
    set.seed(10)
    x_vec <- c(1:n, 1:n)
    y_vec <- c(rna_all[,gene], atac_pred2[,gene])
    col <- c(rep(rgb(0.5, 0.5, 0.5, 0.2), n), rep(rgb(0.8, 0, 0, 0.2), n))
    idx <- sample(1:length(x_vec))
    plot(x = x_vec[idx], y = y_vec[idx], 
         pch = 16, col = col[idx], main = gene)
    lines(rna_all_smoothed[,gene], col = "white", lwd = 4)
    lines(rna_all_smoothed[,gene], col = 1, lwd = 2)
    
    lines(atac_pred2_smoothed[,gene], col = "white", lwd = 4)
    lines(atac_pred2_smoothed[,gene], col = 2, lwd = 2)
  }
  graphics.off()
}

