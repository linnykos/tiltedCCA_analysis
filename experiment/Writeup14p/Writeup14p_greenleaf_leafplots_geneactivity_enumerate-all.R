rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/Writeup14p/Writeup14p_10x_greenleaf_tcca.RData")
load("../../../../out/main/10x_greenleaf_preprocessed_customGAct.RData")

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_common <- multiSVD_obj$common_mat_1[cell_idx[order_vec],]

mentioned_genes <- colnames(rna_common)
# load("../../../../out/main/10x_greenleaf_developmentalGenes.RData")
# mentioned_genes <- c("ALDH2", "APOE", "AQP4", "ASCL1", "BHLHE22",
#                      "BHLHE40", "C16orf89", "CAV2", "DLX2", "DOK5",
#                      "DUSP1", "EOMES", "ETV4", "FOS", "FOXJ1", "GLI3",
#                      "HAS2", "HES1", "HES4", "HSPA1A", "HSPA1B",
#                      "ID3", "IGFBP7", "JUN", "KIF1A", "LIMCH1",
#                      "MBP", "MEF2C", "NEUROD1", "NEUROD2", "NEUROD4",
#                      "NEUROD6", "NEUROG1", "NEUROG2", "NFIA", "NFIB", "NFIC",
#                      "NHLH1", "NR2F1", "PAX6", "RFX4", "RUNX1",
#                      "OLIG1", "OLIG2", "SOX2", "SOX3", "SOX6",
#                      "SOX9", "SOX10", "SOX21", "SPARCL1", "SNCB", "TBX",
#                      "TNC", "TOP2A", "TRB1", "WNT11",
#                      "SATB2")
# mentioned_genes <- unique(c(mentioned_genes, selection_res$selected_variables))
mentioned_genes <- mentioned_genes[which(paste0("ATAC-", mentioned_genes) %in% rownames(greenleaf[["customGAct"]]@scale.data))]
mentioned_genes <- sort(unique(mentioned_genes))
rna_common <- rna_common[,mentioned_genes]

set.seed(10)
svd_atac <- irlba::irlba(greenleaf[["customGAct"]]@scale.data, nv = 50)
atac_mat <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_atac$u[,2:50], svd_atac$d[2:50]), svd_atac$v[,2:50])
rownames(atac_mat) <- rownames(greenleaf[["customGAct"]]@scale.data)
colnames(atac_mat) <- colnames(greenleaf[["customGAct"]]@scale.data)
atac_mat <- atac_mat[paste0("ATAC-",mentioned_genes),cell_idx[order_vec]]
atac_mat <- t(atac_mat)
colnames(atac_mat) <- colnames(rna_common)
all(rownames(atac_mat) == rownames(rna_common))

n <- nrow(rna_common)
p <- ncol(rna_common)
pred_rna_common <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = rna_common[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
pred_atac_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = atac_mat[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(pred_atac_mat) <- colnames(atac_mat)
colnames(pred_rna_common) <- colnames(rna_common)

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))

for(k in 1:ceiling(ncol(rna_common)/30)){
  print(k)
  genes <- mentioned_genes[((k-1)*30+1):min((k*30), length(mentioned_genes))]
  
  png(paste0("../../../../out/figures/Writeup14p/10x_greenleaf_leafplots_geneactivity_enumerate_", 
             k, ".png"),
      height = 3500, width = 3000, units = "px", res = 300)
  par(mfrow = c(6,5), mar = c(4,4,4,0.5))
  for(gene in genes){
    set.seed(10)
    idx <- sample(1:n)
    
    x_vec1 <- atac_mat[,gene]
    x_vec2 <- pred_atac_mat[,gene]
    y_vec1 <- rna_common[,gene]
    y_vec2 <- pred_rna_common[,gene]
    
    x_mean <- mean(x_vec2); x_vec1 <- x_vec1- x_mean; x_vec2 <- x_vec2 - x_mean
    x_sd <- sd(x_vec2); x_vec1 <- x_vec1/x_sd; x_vec2 <- x_vec2/x_sd
    y_mean <- mean(y_vec2); y_vec1 <- y_vec1- y_mean; y_vec2 <- y_vec2 - y_mean
    y_sd <- sd(y_vec2); y_vec1 <- y_vec1/y_sd; y_vec2 <- y_vec2/y_sd
    
    plot(x = x_vec1[idx], y = y_vec1[idx],
         pch = 16, col = color_vec[idx], main = gene, cex = 0.5,
         xlim = quantile(x_vec1, probs = c(0.01, 0.99)),
         ylim = quantile(y_vec1, probs = c(0.01, 0.99)),
         xlab = "ATAC", ylab = "RNA")
    
    points(x = x_vec2[idx], 
           y = y_vec2[idx], 
           pch = 16, col = "white", cex = 2)
    points(x = x_vec2[idx], 
           y = y_vec2[idx], 
           pch = 16, col = color_vec[idx], cex = 1.5)
    
    points(x = x_vec2[1], 
           y = y_vec2[1],
           pch = 16, col = "white", cex = 3)
    points(x = x_vec2[1], 
           y = y_vec2[1], 
           pch = 16, col = color_vec[1], cex = 2.5)
    
    x_right <- mean(x_vec1[x_vec1 >= quantile(x_vec1, probs = 0.9)])
    x_left <- mean(x_vec1[x_vec1 <= quantile(x_vec1, probs = 0.1)])
    y_top <- mean(y_vec1[y_vec1 >= quantile(y_vec1, probs = 0.9)])
    y_bot <- mean(y_vec1[y_vec1 <= quantile(y_vec1, probs = 0.1)])
    lines(c(x_right,x_left), c(y_top, y_bot), lwd = 2, lty = 2)
  }
  graphics.off()
}