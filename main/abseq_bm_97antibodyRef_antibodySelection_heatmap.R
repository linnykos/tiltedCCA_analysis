rm(list=ls())
library(igraph)
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/abseq_bm97Ref_tcca.RData")
load("../../../out/main/abseq_bm97Ref_distinct_differential.RData")
load("../../../out/main/abseq_bm97Ref_varSelect.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

celltype_vec <- factor(bm$ct)
celltype_list <- list("NK cells" = c("CD56brightCD16- NK cells", "CD56dimCD16+ NK cells", "NK T cells"),
                      "CD4" = c("CD4+ cytotoxic T cells", "CD4+ memory T cells", "CD4+ naive T cells",
                                "CD69+PD-1+ memory CD4+ T cells"),
                      "CD8" = c("CD8+ central memory T cells", "CD8+ effector memory T cells", "CD8+ naive T cells",
                                "CD8+CD103+ tissue resident memory T cells"), 
                      "B" = c("CD11c+ memory B cells", "Class switched memory B cells", "Immature B cells",
                              "Mature naive B cells", "Nonswitched memory B cells", "Pre-B cells",
                              "Pro-B cells"),
                      "Erythroid" = c("Early erythroid progenitor",
                                      "Eosinophil-basophil-mast cell progenitors", "Erythro-myeloid progenitors",
                                      "Late erythroid progenitor", "Megakaryocyte progenitors"),
                      "HSC-MPPs" = c("HSCs & MPPs"),
                      "Myeloid" = c("Classical Monocytes", "Conventional dendritic cell 1",
                                    "Conventional dendritic cell 1", "Early promyelocytes",
                                    "Late promyelocytes", "Lymphomyeloid prog", "Myelocytes",
                                    "Non-classical monocytes", "Plasmacytoid dendritic cell progenitors",
                                    "Plasmacytoid dendritic cells"),
                      "Other" = c("GammaDelta T cells", "Plasma cells", "Aberrant erythroid"))
celltype_all <- unlist(celltype_list)
names(celltype_all) <- NULL

protein_mat <- as.matrix(Matrix::t(bm[["AB"]]@scale.data))
protein_mat <- protein_mat[,variable_selection_res$selected_variables]
protein_avg <- t(sapply(celltype_all, function(celltype){
  idx <- which(celltype_vec == celltype)
  colMeans(protein_mat[idx,,drop = F])
}))
protein_avg <- scale(protein_avg)

# hclust_res <- stats::hclust(stats::dist(t(protein_avg)))
# protein_avg2 <- protein_avg[,hclust_res$order]

break_vec <- c(seq(min(protein_avg), 0, length.out = 11), seq(0, max(protein_avg), length.out = 11)[-1])
break_vec[1] <- break_vec[1]-1
break_vec[length(break_vec)] <- break_vec[length(break_vec)]+1
col_vec1 <- grDevices::colorRampPalette(c(rgb(184, 54, 220, maxColorValue = 255), "white"))(11)[-1]
col_vec2 <- grDevices::colorRampPalette(c("white",  rgb(235, 134, 47, maxColorValue = 255)))(11)[-1]
col_vec <- c(col_vec1, col_vec2)

png("../../../out/figures/main/abseq_bm97Ref_varSelect_protein-everything_heatmap.png",
    height = 3500, width = 3500*ncol(protein_avg)/nrow(protein_avg), res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(protein_avg), 
      asp = (nrow(protein_avg)-1)/(ncol(protein_avg)-1),
      main = "", xlab = "", ylab = "",
      xaxt = "n", yaxt = "n", bty = "n",
      breaks = break_vec, col = col_vec)

# draw lines
n <- nrow(protein_avg)
halfspacing <- 1/(2*(n-1))
num_types <- sapply(celltype_list, length)
for(i in 1:(length(num_types)-1)){
  y_val <- 1-((2*halfspacing)*(sum(num_types[1:i])-1)+halfspacing)
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 5, col = "white")
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 4, lty = 2, col = "black")
}
graphics.off()
