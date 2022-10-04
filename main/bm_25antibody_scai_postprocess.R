rm(list=ls())
library(Seurat)

load("../../../out/main/citeseq_bm25_scAI.RData")
source("bm_25antibody_colorPalette.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

########

set.seed(10)
umap_res <- Seurat::RunUMAP(scai_res$H)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(bm)
colnames(umap_mat) <- paste0("scaiUMAP_", 1:ncol(umap_mat))
bm[["scaiUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                 assay = "RNA")

plot1 <- Seurat::DimPlot(bm, reduction = "scaiUMAP",
                         group.by = "celltype.l2", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_bm25_scAI_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

#######################


celltype_vec <- factor(bm$celltype.l2)
celltype_list <- list("CD4" = c("CD4 Memory", "CD4 Naive", "Treg"),
                      "CD8" = c("CD8 Effector_1", "CD8 Effector_2", "CD8 Memory_1", "CD8 Memory_2", "CD8 Naive", "gdT", "MAIT"), 
                      "NK cells" = c("CD56 bright NK", "NK"),
                      "B" = c("Memory B", "Naive B", "Plasmablast"),
                      "Myeloid" = c("CD14 Mono", "CD16 Mono", "cDC2", "pDC"),
                      "Progenitor" = c("GMP", "HSC", "LMPP", "Prog_DC", "Prog_Mk", "Prog_RBC", "Prog_B 1", "Prog_B 2"))
celltype_all <- unlist(celltype_list)
names(celltype_all) <- NULL


.processing_func <- function(mat, val = 0.05){
  mat <- scale(mat, center = T, scale = F)
  pmin(pmax(mat, quantile(mat, probs = val)), quantile(mat,probs = 1-val))
}

scai_mat <- scai_res$H
scai_mat <- .processing_func(scai_mat)
scai_avg <- t(sapply(celltype_all, function(celltype){
  idx <- which(celltype_vec == celltype)
  colMeans(scai_mat[idx,,drop = F])
}))

break_vec <- c(seq(min(scai_avg), 0, length.out = 11), seq(0, max(scai_avg), length.out = 11)[-1])
break_vec[1] <- break_vec[1]-1
break_vec[length(break_vec)] <- break_vec[length(break_vec)]+1
base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
col_vec <- rev(grDevices::colorRampPalette(base_palette)(length(break_vec)-1))

png("../../../out/figures/main/bm_25antibody_scai-heatmap.png",
    height = 3500, width = 3500*ncol(scai_avg)/nrow(scai_avg), res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(scai_avg), 
      asp = (nrow(scai_avg)-1)/(ncol(scai_avg)-1),
      main = "", xlab = "", ylab = "",
      xaxt = "n", yaxt = "n", bty = "n",
      breaks = break_vec, col = col_vec)

# draw lines
n <- nrow(scai_avg)
halfspacing <- 1/(2*(n-1))
num_types <- sapply(celltype_list, length)
for(i in 1:(length(num_types)-1)){
  y_val <- 1-((2*halfspacing)*(sum(num_types[1:i])-1)+halfspacing)
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 5, col = "white")
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 4, lty = 2, col = "black")
}

# draw lines
n <- ncol(scai_avg)
halfspacing <- 1/(2*(n-1))
num_types <- rep(5,10)
for(j in 1:(length(num_types)-1)){
  y_val <- 1-((2*halfspacing)*(sum(num_types[1:j])-1)+halfspacing)
  lines(y = c(0,1), x = rep(y_val, 2), lwd = 5, col = "white")
  lines(y = c(0,1), x = rep(y_val, 2), lwd = 4, lty = 2, col = "black")
}

graphics.off()

##############

rna_loadings <- apply(scai_res$W1, 2, mean)/quantile(scai_res$W1, prob = 0.95)
adt_loadings <- apply(scai_res$W2, 2, mean)/quantile(scai_res$W2, prob = 0.95)

df <- rbind(rna_loadings, adt_loadings)

png("../../../out/figures/main/bm_25antibody_scai_factor-weights.png",
    height = 1000, width = 4000, res = 500, units = "px")
par(mar = c(0.5, 4, 0.5, 0.5))
barplot(df, main = "", ylab = "", xlab = "",
        col = c(rgb(55,100,234, maxColorValue = 255),
                rgb(82,185,44, maxColorValue = 255)),
        beside=TRUE)
graphics.off()


celltype_vec <- bm$celltype.l2
col_vec <- sapply(celltype_vec, function(celltype){
  col_palette[which(names(col_palette) == celltype)]
})

dim_vec <- c(4, 16, 19, 20)
set.seed(10)
idx <- sample(1:nrow(scai_res$H))
for(j in dim_vec){
  vec1 <- scai_res$H[,j]
  png(paste0("../../../out/figures/main/bm_25antibody_scai_factor-",j ,".png"),
      height = 1200, width = 1000, units = "px", res = 500)
  par(mar = c(0.5, 2, 0.5, 0.5), bg = NA)
  plot(NA, ylim = range(vec1), xlim = c(-.5, 0.5), 
       ylab = "", xlab = "", bty = "n", xaxt = "n")
  axis(side = 2, labels = F)
  points(y = vec1[idx], x = rep(0, length(vec1)) + runif(length(vec1), min = -.3, max = .3),
         pch = 16, col = col_vec[idx], cex = 0.25)
  graphics.off()

}

