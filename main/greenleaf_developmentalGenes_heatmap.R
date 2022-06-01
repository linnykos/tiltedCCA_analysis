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
for(i in 1:nrow(variable_summary)){
  if(variable_summary[i,"pseudotime"] <= 0.75 & variable_summary[i,"Lineage"] == 2){
    variable_summary[i,"Lineage"] <- 1
  }
}

# extract the relevant matrix
lineage_idx_all <- sort(unique(c(lineage_idx1, lineage_idx2)))
heatmap_mat <- greenleaf[["SCT"]]@data[selected_variables,lineage_idx_all]
heatmap_mat <- as.matrix(heatmap_mat)
heatmap_mat <- scale(t(heatmap_mat))

# normalize the cells in the matrix

# order the genes in the matrix

# split the matrix by cells per-lineage