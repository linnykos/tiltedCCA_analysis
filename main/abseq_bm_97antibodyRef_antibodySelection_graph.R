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

Seurat::DefaultAssay(bm) <- "AB"
mat_2 <- Matrix::t(bm[["AB"]]@scale.data)
mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}
cor_mat <- abs(stats::cor(mat_2b))
diag(cor_mat) <- 0

# compute the layout via a graph
n <- ncol(cor_mat)
adj_mat <- matrix(0, n, n)
# adj_mat[cor_mat >= variable_selection_res$cor_threshold] <- cor_mat[cor_mat >= variable_selection_res$cor_threshold]
adj_mat[cor_mat >= 0.7] <- 1
k <- 10
for(i in 1:n){
  idx <- order(cor_mat[,i], decreasing = T)[1:k]
  val_vec <- sort(cor_mat[,i], decreasing = T)[1:k]
  adj_mat[i,idx] <- 1
}
adj_mat <- (adj_mat + t(adj_mat))/2
rownames(adj_mat) <- colnames(mat_2b)
colnames(adj_mat) <- colnames(mat_2b)
set.seed(10)
umap_res <- Seurat::RunUMAP(adj_mat, spread = 5, min.dist = 0.1)

cor_power <- 5; max_width <- 3
g <- igraph::graph.adjacency(adj_mat, mode="undirected")
# igraph::E(g)$width <- igraph::E(g)$weight^cor_power
# igraph::E(g)$width <- max_width*igraph::E(g)$width/max(igraph::E(g)$width)

g2 <- g
color_palette <- grDevices::colorRampPalette(c("white", 2))(20)
cor_power <- 1; max_val <- 100
x_vec <- logpval_vec; x_vec[x_vec >= max_val] <- 2*max_val
x_vec <- (x_vec - min(x_vec))/diff(range(x_vec)); x_vec <- x_vec^cor_power
seq_vec <- seq(0, 1, length.out = 20)
col_vec <- sapply(names(igraph::V(g2)), function(i){
  val <- x_vec[i]
  idx <- which.min(abs(seq_vec - val))
  color_palette[idx]
})
igraph::V(g2)$color <- col_vec

size_breakpoints <- c(min(rsquare_vec)-.1, 0.75, 0.84, max(rsquare_vec)+.1)
size_vals <- c(8,6,4)
size_vec <- sapply(names(igraph::V(g2)), function(i){
  val <- rsquare_vec[i]
  idx <- max(which(size_breakpoints <= val))
  size_vals[idx]
})
igraph::V(g2)$size <- size_vec

png("../../../out/figures/main/abseq_bm97Ref_varSelect_graph-cleaned.png",
    height = 1500, width = 2000, units = "px", res = 500)
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

png("../../../out/figures/main/abseq_bm97Ref_varSelect_graph-labeled.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = rep(0.5,4))
set.seed(10)
plot(g2, 
     layout = umap_res@cell.embeddings,
     vertex.size = 3, 
     edge.curved = 0.3,
     vertex.label.cex = 0.5)
graphics.off()

###############################

variable_selection_res$selected_variables
table(adj_mat["Disialoganglioside.GD2-AB",])
vec <- adj_mat["CD40-AB",]; vec[which(vec != 0)]
vec <- adj_mat["IgD-AB",]; vec[which(vec != 0)]

################################3

multiSVD_obj2 <- tiltedCCA:::tiltedCCA_decomposition(multiSVD_obj, 
                                     bool_modality_1_full = F,
                                     bool_modality_2_full = F,
                                     verbose = verbose)
reference_dimred <- tiltedCCA:::.get_tCCAobj(multiSVD_obj2, apply_postDimred = F, what = "common_dimred")
# reference_dimred <- cbind(reference_dimred)#, mat_2b[,variable_selection_res$selected_variables[1]])
tiltedCCA:::.linear_regression(bool_include_intercept = T,
                   bool_center_x = T,
                   bool_center_y = T,
                   bool_scale_x = T,
                   bool_scale_y = T,
                   return_type = "r_squared", 
                   x_mat = reference_dimred,
                   y_vec = mat_2b[,"CD8-AB"])