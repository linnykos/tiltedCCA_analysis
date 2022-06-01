rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")
load("../../../out/main/10x_greenleaf_developmentalGenes.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

selected_variables <- selection_res$selected_variables

# for each gene, fit two np regression based on Lineage1 and Lineage2, and assign 
#  the pseudotime and branch for each gene
lineage_idx1 <- which(greenleaf$Lineage1 == 1)
lineage_idx2 <- which(greenleaf$Lineage2 == 1)
ordering1 <- order(greenleaf$pseudotime[lineage_idx1], decreasing = F)
ordering2 <- order(greenleaf$pseudotime[lineage_idx2], decreasing = F)
x_vec1 <- greenleaf$pseudotime[lineage_idx1[ordering1]]
x_vec2 <- greenleaf$pseudotime[lineage_idx2[ordering2]]
names(x_vec1) <- NULL
names(x_vec2) <- NULL

variable_summary <- sapply(1:length(selected_variables), function(i){
  print(paste0("Working on gene ", i, " out of ", length(selected_variables)))
  gene <- selected_variables[i]
  
  y_vec1 <- greenleaf[["SCT"]]@data[gene,lineage_idx1[ordering1]]
  y_vec2 <- greenleaf[["SCT"]]@data[gene,lineage_idx2[ordering2]]
  names(y_vec1) <- NULL
  names(y_vec2) <- NULL
  
  df1 <- data.frame(x = x_vec1,  y = y_vec1)
  np_res1 <- npregfast::frfast(y ~ x, data = df1)
  max_val1 <- max(np_res1$p[,1,1])
  pseudotime1 <- np_res1$x[which.max(np_res1$p[,1,1])]
  
  df2 <- data.frame(x = x_vec2,  y = y_vec2)
  np_res2 <- npregfast::frfast(y ~ x, data = df2)
  max_val2 <- max(np_res2$p[,1,1])
  pseudotime2 <- np_res2$x[which.max(np_res2$p[,1,1])]
  
  if(max_val1 > max_val2){
    return(c(Lineage = 1, pseudotime = pseudotime1))
  } else {
    return(c(Lineage = 2, pseudotime = pseudotime2))
  }
})
variable_summary <- t(variable_summary)
rownames(variable_summary) <- selected_variables

sapply(unique(greenleaf$celltype), function(celltype){
  vec <- greenleaf$pseudotime[which(greenleaf$celltype == celltype)]
  vec <- vec[!is.na(vec)]
  stats::quantile(vec)
})

# hard set pseudotimes too early as lineage 1 for visualization purposes
branching_threshold <- 0.75
for(i in 1:nrow(variable_summary)){
  if(variable_summary[i,"pseudotime"] <= branching_threshold & variable_summary[i,"Lineage"] == 2){
    variable_summary[i,"Lineage"] <- 1
  }
}

# extract the relevant matrix
lineage_idx_all <- sort(unique(c(lineage_idx1, lineage_idx2)))
heatmap_mat <- greenleaf[["SCT"]]@data[selected_variables,lineage_idx_all]
heatmap_mat <- as.matrix(heatmap_mat)

# normalize the cells in the matrix
heatmap_mat <- scale(t(heatmap_mat))

# order the genes in the matrix
tmp1 <- which(variable_summary[,"Lineage"] == 1)
gene_ordering1 <- rownames(variable_summary)[tmp1[order(variable_summary[tmp1,"pseudotime"], decreasing = F)]]
tmp2 <- which(variable_summary[,"Lineage"] == 2)
gene_ordering2 <- rownames(variable_summary)[tmp2[order(variable_summary[tmp2,"pseudotime"], decreasing = F)]]
gene_ordering <- c(gene_ordering1, gene_ordering2)
heatmap_mat <- heatmap_mat[,gene_ordering]

# split the matrix by cells per-lineage
cell_names <- rownames(heatmap_mat)
tmp <- which(colnames(greenleaf) %in% cell_names)
cell_names2 <- colnames(greenleaf)[tmp]
df_mat <- t(sapply(1:length(tmp), function(i){
  if(i %% floor(length(tmp)/10) == 0) cat('*')
  
  idx <- tmp[i]
  pseudotime <- greenleaf$pseudotime[idx]; names(pseudotime) <- NULL
  bool_lineage1 <- (greenleaf$Lineage1[idx] == 1)
  if(pseudotime < branching_threshold | bool_lineage1){
    return(c(Lineage = 1, pseudotime = pseudotime))
  } else {
    return(c(Lineage = 2, pseudotime = pseudotime))
  }
}))
rownames(df_mat) <- cell_names2
tmp1 <- which(df_mat[,"Lineage"] == 1)
cell_ordering1 <- rownames(df_mat)[tmp1[order(df_mat[tmp1,"pseudotime"], decreasing = F)]]
tmp2 <- which(df_mat[,"Lineage"] == 2)
cell_ordering2 <- rownames(df_mat)[tmp2[order(df_mat[tmp2,"pseudotime"], decreasing = F)]]
cell_ordering <- c(cell_ordering1, cell_ordering2)


# set up the color palette
heatmap_mat1 <- t(heatmap_mat[cell_ordering1,])
heatmap_mat1 <- sign(heatmap_mat1) * abs(heatmap_mat1)^0.5 #^scaling_grid[which.max(scaling_quality)]
heatmap_mat2 <- t(heatmap_mat[cell_ordering2,])
heatmap_mat2 <- sign(heatmap_mat2) * abs(heatmap_mat2)^0.5 # ^scaling_grid[which.max(scaling_quality)]

# remove ROBO2
heatmap_mat1 <- heatmap_mat1[-which(rownames(heatmap_mat1) == "ROBO2"),]
heatmap_mat2 <- heatmap_mat2[-which(rownames(heatmap_mat2) == "ROBO2"),]

num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(c(heatmap_mat1, heatmap_mat2)), 
                    max(c(heatmap_mat1, heatmap_mat2)), 
                    length.out = num_color+1)

png(paste0("../../../out/figures/main/10x_greenleaf_developmentalGenes_heatmap1.png"),
    height = 1300, width = 2000, units = "px", res = 500)
par(mar = c(0,0,0,0))
image(tiltedCCA:::.rotate(heatmap_mat1), xaxs = "i", yaxs = "i",
      col = color_palette, breaks = color_breaks,
      main = "",
      xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")

y_val <- 1-length(gene_ordering1)/length(gene_ordering)
lines(x = c(0,1), y = rep(y_val, 2), col = "white", lwd = 3)
lines(x = c(0,1), y = rep(y_val, 2), lty = 2, lwd = 2)
graphics.off()

png(paste0("../../../out/figures/main/10x_greenleaf_developmentalGenes_heatmap2.png"),
    height = 1300, width = ncol(heatmap_mat2)*2000/ncol(heatmap_mat1), units = "px", res = 500)
par(mar = c(0,0,0,0))
image(tiltedCCA:::.rotate(heatmap_mat2), xaxs = "i", yaxs = "i",
      col = color_palette, breaks = color_breaks,
      main = "",
      xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")

y_val <- 1-length(gene_ordering1)/length(gene_ordering)
lines(x = c(0,1), y = rep(y_val, 2), col = "white", lwd = 3)
lines(x = c(0,1), y = rep(y_val, 2), lty = 2, lwd = 2)
graphics.off()

##################################

# steady state plot
rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(rna = rna_common[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})
greenleaf$alignment <- alignment_vec

lineage_idx1 <- which(greenleaf$Lineage1 == 1)
lineage_idx2 <- which(greenleaf$Lineage2 == 1)
