rm(list=ls())
load("../../../out/main/citeseq_pbmc224_differential.RData")
load("../../../out/main/citeseq_pbmc224_varSelect.RData")
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")
library(igraph)

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

############################################

Seurat::DefaultAssay(pbmc) <- "ADT"
logpval_vec <- sapply(rownames(pbmc), function(gene){
  val <- stats::quantile(
    sapply(1:length(adt_de_list), function(j){
      idx <- which(rownames(adt_de_list[[j]]) == gene)
      if(length(idx) == 0) return(1)
      adt_de_list[[j]][idx, "p_val"]
    }), probs = 0.75
  )
  
  -log10(val)
})
# tmp <- sapply(1:length(adt_de_list), function(j){
#   adt_de_list[[j]][which(rownames(adt_de_list[[j]]) == "CD52"), "p_val"]
# })
names(logpval_vec) <- rownames(pbmc)
quantile(logpval_vec)
max_val <- max(logpval_vec[which(!is.infinite(logpval_vec))])
logpval_vec <- pmin(logpval_vec, 300)

#######################

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(multiSVD_obj, 
                                                    bool_modality_1_full = F,
                                                    bool_modality_2_full = F,
                                                    verbose = 0)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
reference_dimred <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "common_dimred")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
distinct_mat <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "distinct_mat")

rsquare_vec <- sapply(1:ncol(distinct_mat), function(j){
  print(j)
  tiltedCCA:::.linear_regression(bool_include_intercept = T,
                                 bool_center_x = T,
                                 bool_center_y = T,
                                 bool_scale_x = T,
                                 bool_scale_y = T,
                                 return_type = "r_squared", 
                                 x_mat = reference_dimred,
                                 y_vec = distinct_mat[,j])
})
names(rsquare_vec) <- colnames(distinct_mat)
rsquare_vec <- rsquare_vec[colnames(multiSVD_obj$distinct_mat_2)]

############################################

cor_power <- 4

mat <- multiSVD_obj$distinct_mat_2
cor_mat <- stats::cor(mat)
cor_mat2 <- abs(cor_mat)^cor_power
diag(cor_mat2) <- 0

g <- igraph::graph.adjacency(cor_mat2, mode="undirected", weighted=TRUE)
igraph::E(g)$width <- igraph::E(g)$weight*4

g2 <- g
x_vec <- rsquare_vec
color_palette <- grDevices::colorRampPalette(c("white", 2))(20)
seq_vec <- seq(min(x_vec), max(x_vec), length.out = 20)
col_vec <- sapply(names(igraph::V(g2)), function(i){
  val <- x_vec[i]
  idx <- which.min(abs(seq_vec - val))
  color_palette[idx]
})
igraph::V(g2)$color <- col_vec
png("../../../out/figures/main/citeseq_pbmc224_varSelect_graph-rsquared.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = rep(0.5,4))
set.seed(10)
plot(g2, 
     vertex.size = 3, 
     vertex.label = "", 
     vertex.color = igraph::V(g2)$color, 
     edge.curved = 0.3)
graphics.off()

############

g2 <- g
x_vec <- sqrt(logpval_vec)
color_palette <- grDevices::colorRampPalette(c("white", 2))(20)
seq_vec <- seq(min(x_vec), max(x_vec), length.out = 20)
col_vec <- sapply(names(igraph::V(g2)), function(i){
  val <- x_vec[i]
  idx <- which.min(abs(seq_vec - val))
  color_palette[idx]
})
igraph::V(g2)$color <- col_vec
png("../../../out/figures/main/citeseq_pbmc224_varSelect_graph-logpval.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = rep(0.5,4))
set.seed(10)
plot(g2, 
     vertex.size = 3, 
     vertex.label = "", 
     vertex.color = igraph::V(g2)$color, 
     edge.curved = 0.3)
graphics.off()

############

g2 <- g
col_vec <- sapply(names(igraph::V(g2)), function(i){
  if(i %in% variable_selection_res$selected_variables) "deepskyblue2" else "white"
})
igraph::V(g2)$color <- col_vec
png("../../../out/figures/main/citeseq_pbmc224_varSelect_graph-selected.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = rep(0.5,4))
set.seed(10)
plot(g2, 
     vertex.size = 3, 
     vertex.label = "", 
     vertex.color = igraph::V(g2)$color, 
     edge.curved = 0.3)
graphics.off()


############

g2 <- g
igraph::V(g2)$label <- names(rsquare_vec)
png("../../../out/figures/main/citeseq_pbmc224_varSelect_graph-labeled.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = rep(0.5,4))
set.seed(10)
plot(g2, 
     vertex.size = 3, 
     vertex.cex = 0.3,
     vertex.label = igraph::V(g2)$label, 
     edge.curved = 0.3)
graphics.off()


############

g2 <- g
igraph::V(g2)$color <- col_vec
png("../../../out/figures/main/citeseq_pbmc224_varSelect_graph-selected.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = rep(0.5,4))
set.seed(10)
plot(g2, 
     vertex.size = 3, 
     vertex.label = "", 
     vertex.color = igraph::V(g2)$color, 
     edge.curved = 0.3)
graphics.off()

###############

logpval_vec["CD18"]
rsquare_vec["CD18"]

logpval_vec[variable_selection_res$selected_variables]
rsquare_vec[variable_selection_res$selected_variables]

for(i in 1:length(variable_selection_res$candidate_list)){
  print("=====")
  print(i)
  print(sort(variable_selection_res$candidate_list[[i]]))
}