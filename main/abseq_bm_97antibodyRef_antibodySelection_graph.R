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

##############3

Seurat::DefaultAssay(bm) <- "AB"
antibody_names <- rownames(bm)
logpval_vec <- sapply(1:length(antibody_names), function(k){
  if(k %% floor(length(antibody_names)/10) == 0) cat('*')
  antibody <- antibody_names[k]
  
  # cycle through all the celltypes
  celltype_vec <- sapply(1:length(adt_distinct_de_list$level_vec), function(i){
    idx <- which(adt_distinct_de_list$combn_mat == i, arr.ind = T)[,2]
    vec <-  sapply(idx, function(j){
      idx <- which(rownames(adt_distinct_de_list$de_list[[j]]) == antibody)
      if(length(idx) == 0) return(1)
      adt_distinct_de_list$de_list[[j]][idx, "p_val"]
    })
    stats::quantile(vec, probs = 0.75)
  })
  
  max(-log10(celltype_vec))
})
logpval_vec <- pmin(logpval_vec, 300)
names(logpval_vec) <- rownames(bm)

rsquare_vec <- variable_selection_res$cor_vec_intial
logpval_vec <- logpval_vec[names(logpval_vec) %in% names(rsquare_vec)]
rsquare_vec <- rsquare_vec[names(logpval_vec)]

logpval_vec[variable_selection_res$selected_variables]
rsquare_vec[variable_selection_res$selected_variables]

quantile(logpval_vec)
quantile(rsquare_vec)

##########################

cor_power <- 3
Seurat::DefaultAssay(bm) <- "AB"
mat_2 <- Matrix::t(bm[["AB"]]@scale.data)
mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}
cor_mat <- abs(stats::cor(mat_2b))
diag(cor_mat) <- 0
cor_mat2 <- cor_mat^cor_power

set.seed(10)
umap_res <- Seurat::RunUMAP(cor_mat2)

cor_power <- 5; max_width <- 3
g <- igraph::graph.adjacency(cor_mat, mode="undirected", weighted=TRUE)
igraph::E(g)$width <- igraph::E(g)$weight^cor_power
igraph::E(g)$width <- max_width*igraph::E(g)$width/max(igraph::E(g)$width)

g2 <- g
color_palette <- grDevices::colorRampPalette(c("white", 2))(20)
cor_power <- 1
x_vec <- logpval_vec; x_vec <- (x_vec - min(x_vec))/diff(range(x_vec)); x_vec <- x_vec^cor_power
seq_vec <- seq(0, 1, length.out = 20)
col_vec <- sapply(names(igraph::V(g2)), function(i){
  val <- x_vec[i]
  idx <- which.min(abs(seq_vec - val))
  color_palette[idx]
})
igraph::V(g2)$color <- col_vec

size_breakpoints <- c(min(rsquare_vec)-.1, 0.8, 0.85, max(rsquare_vec)+.1)
size_vals <- c(4,2,1)
size_vec <- sapply(names(igraph::V(g2)), function(i){
  val <- rsquare_vec[i]
  idx <- max(which(size_breakpoints <= val))
  size_vals[idx]
})
igraph::V(g2)$size <- size_vec

png("../../../out/figures/main/abseq_bm97Ref_varSelect_graph-cleaned.png",
    height = 3000, width = 3000, units = "px", res = 500)
par(mar = rep(0.5,4))
set.seed(10)
plot(g2, 
     layout = umap_res@cell.embeddings,
     vertex.size = igraph::V(g2)$size, 
     vertex.label = "", 
     vertex.color = igraph::V(g2)$color, 
     edge.curved = 0.3)
graphics.off()

####

g2 <- g
col_vec <- sapply(names(igraph::V(g2)), function(i){
  if(i %in% variable_selection_res$selected_variables) "deepskyblue2" else "white"
})
igraph::V(g2)$color <- col_vec
png("../../../out/figures/main/abseq_bm97Ref_varSelect_graph-selected.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = rep(0.5,4))
set.seed(10)
plot(g2, 
     layout = umap_res@cell.embeddings,
     vertex.size = 3, 
     vertex.label = "", 
     vertex.color = igraph::V(g2)$color, 
     edge.curved = 0.3)
graphics.off()
