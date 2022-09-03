rm(list=ls())
library(Seurat); library(Signac)
load("../../../out/main/10x_greenleaf_tcca_RNA-geneActivity.RData")
source("greenleaf_colorPalette.R")

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+Gene Activity)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-geneActivity_umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::FeaturePlot(greenleaf, reduction = "common_tcca",
                             features = "pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+Gene Activity): T-CCA's Common\nSlingshot's pseudotime via ATAC"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-geneActivity_umap_common-pseudotime.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(greenleaf, reduction = "distinct1_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+Gene Activity)\nRNA distinct subspace"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-geneActivity_umap_distinct1.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(greenleaf, reduction = "distinct2_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+Gene Activity)\nGene Activity distinct subspace"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-geneActivity_umap_distinct2.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

##########################

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "celltype", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_geneActivity-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

set.seed(10)
greenleaf <- Seurat::FindMultiModalNeighbors(greenleaf, 
                                             reduction.list = list("pca", "pcaCustomGAct"), 
                                             weighted.nn.name = "weighted.nn2",
                                             dims.list = list(1:50, 2:50))
set.seed(10)
greenleaf <- Seurat::RunUMAP(greenleaf, nn.name = "weighted.nn2", reduction.name = "umap.wnn2", 
                             reduction.key = "wnnUMAP2_")

plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap.wnn2",
                         group.by = "celltype", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_wnn-umap_RNA-geneActivity_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

# consensus pca
Seurat::DefaultAssay(greenleaf) <- "SCT"
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[Seurat::VariableFeatures(object = greenleaf),])
Seurat::DefaultAssay(greenleaf) <- "customGAct"
mat_2 <- Matrix::t(greenleaf[["customGAct"]]@data[Seurat::VariableFeatures(object = greenleaf),])

mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = mat_1b, mat_2 = mat_2b,
                                           dims_1 = 1:50, dims_2 = 2:50,
                                           dims_consensus = 1:49,
                                           center_1 = T, center_2 = T,
                                           recenter_1 = F, recenter_2 = F,
                                           rescale_1 = F, rescale_2 = F,
                                           scale_1 = T, scale_2 = T,
                                           verbose = 1)
consensus_dimred <- consensus_pca$dimred_consensus
colnames(consensus_dimred) <- paste0("consensusPCA2_", 1:ncol(consensus_dimred))
greenleaf[["consensusPCA2"]] <- Seurat::CreateDimReducObject(consensus_dimred, 
                                                             assay = "SCT")

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(greenleaf)
colnames(umap_mat) <- paste0("consensusUMAP2_", 1:ncol(umap_mat))
greenleaf[["consensusUMAP2"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                              assay = "SCT")

plot1 <- Seurat::DimPlot(greenleaf, reduction = "consensusUMAP2",
                         group.by = "celltype", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_consensusPCA-umap_RNA-geneActivity_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

##########################

var_features <- sort(c("ALDH2", "APOE", "AQP4", "ASCL1", "BHLHE22",
                       "BHLHE40", "C16orf89", "CAV2", "DLX2", "DOK5",
                       "DUSP1", "EOMES", "ETV4", "FOS", "FOXJ1", "GLI3",
                       "HAS2", "HES1", "HES4", "HSPA1A", "HSPA1B",
                       "ID3", "IGFBP7", "JUN", "KIF1A", "LIMCH1",
                       "MBP", "MEF2C", "NEUROD1", "NEUROD2", "NEUROD4",
                       "NEUROD6", "NEUROG1", "NEUROG2", "NFIA", "NFIB", "NFIC",
                       "NHLH1", "NR2F1", "PAX6", "RFX4", "RUNX1",
                       "OLIG1", "OLIG2", "SOX2", "SOX3", "SOX6",
                       "SOX9", "SOX10", "SOX21", "SPARCL1", "SNCB", "TBX",
                       "TNC", "TOP2A", "TRB1", "WNT11",
                       "SOX4", "DSCAM", "SATB2", "DAB1", "OPCML", "GRIN2B",
                       "MDGA2", "MEIS2", "CNTN4", "TENM2", "ADGRL3",
                       "VSTM2L", "CHL1"))
var_features <- intersect(var_features, rownames(greenleaf[["SCT"]]))
var_features <- var_features[which(paste0("ATAC-",var_features) %in% rownames(greenleaf[["customGAct"]]))]

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_mat <- multiSVD_obj$common_mat_1[cell_idx[order_vec],var_features]
atac_mat <- multiSVD_obj$common_mat_2[cell_idx[order_vec],paste0("ATAC-", var_features)]

n <- nrow(rna_mat)
p <- ncol(rna_mat)
pred_rna_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = rna_mat[,j], x = 1:n)
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
colnames(pred_rna_mat) <- colnames(rna_mat)

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))

for(k in 1:ceiling(ncol(rna_mat)/20)){
  print(k)
  genes <- colnames(rna_mat)[((k-1)*20+1):min((k*20), ncol(rna_mat))]
  
  png(paste0("../../../out/figures/main/10x_greenleaf_leafplots_RNA-geneActivity_common_enumerate_", 
             k, ".png"),
      height = 3500, width = 3000, units = "px", res = 300)
  par(mfrow = c(5,4), mar = c(4,4,4,0.5))
  for(gene in genes){
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

#########################

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

n <- nrow(rna_mat)
p <- ncol(rna_mat)
pred_rna_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = rna_mat[,j], x = 1:n)
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
colnames(pred_rna_mat) <- colnames(rna_mat)

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))

for(k in 1:ceiling(ncol(rna_mat)/20)){
  print(k)
  genes <- colnames(rna_mat)[((k-1)*20+1):min((k*20), ncol(rna_mat))]
  
  png(paste0("../../../out/figures/main/10x_greenleaf_leafplots_RNA-geneActivity_original_enumerate_", 
             k, ".png"),
      height = 3500, width = 3000, units = "px", res = 300)
  par(mfrow = c(5,4), mar = c(4,4,4,0.5))
  for(gene in genes){
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

