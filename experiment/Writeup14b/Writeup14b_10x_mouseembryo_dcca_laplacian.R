rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
source("../Writeup14/Writeup14_peakcalling_function.R")
date_of_run <- Sys.time(); session_info <- sessionInfo()

mat_1 <- t(mbrain[["SCT"]]@scale.data)
Seurat::DefaultAssay(mbrain) <- "ATAC"
mat_2 <- Matrix::t(mbrain[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain),])

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
metadata <- mbrain@meta.data

#######

set.seed(10)
mbrain <- FindMultiModalNeighbors(mbrain, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
mbrain <- FindClusters(mbrain, graph.name = "wsnn", algorithm = 3, resolution=22,verbose = FALSE)
ext.upstream=500000
ext.downstream=500000
# include some cell features object, to aggregate over.
include.cell.features=c("nCount_RNA", "nCount_ATAC", "wsnn_res.22")
gene_vec <- c("Pax6", "Neurog2", "Eomes", "Neurod1", "Neurod2", "Neurod4", "Neurod6", "Tbr1",
              "Sox5", "Fezf2", "Ikzf1", "Foxg1", "Tle4", "Bcl11b", "Nr2f1", "Tbr1",
              "Satb2", "Pou3f2", "Pou3f3")
gene_vec <- gene_vec[gene_vec %in% colnames(mat_1)]
myobj <- getPeaksForGenes(mbrain, gene.names = gene_vec, 
                          include.cell.features=include.cell.features,
                          ext.upstream=ext.upstream, ext.downstream=ext.downstream)
myobj <- aggregateCells(myobj, aggregate.over="wsnn_res.22")
myobj <- findLinks(myobj)

#######

cell_idx <- which(metadata$label_Savercat %in% c("Oligodendrocyte", "Hindbrain glycinergic", "Midbrain glutamatergic",
                                                 "Forebrain glutamatergic", "Radial glia", "Forebrain GABAergic", 
                                                 "Neuroblast", "Cajal-Retzius", "Mixed region GABAergic", 
                                                 "Glioblast", "Cortical or hippocampal glutamatergic"))
set.seed(10)
rank_1 <- 30; rank_2 <- 31
mat_1 <- mat_1[cell_idx,]; mat_2 <- mat_2[cell_idx,]
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = T,
                                      meta_clustering = NA, num_neigh = 15, 
                                      apply_shrinkage = F, fix_distinct_perc = F, 
                                      verbose = T) 
dcca_res2 <- dcca_res
class(dcca_res2) <- "dcca_decomp"

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)
mat_1_denoised <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
mat_2_denoised <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2

############

mbrain2 <- Seurat::CreateSeuratObject(counts = t(mat_1_denoised))
mbrain2[["label_Savercat"]] <- metadata$label_Savercat[cell_idx]
membership_vec <- as.factor(metadata$label_Savercat[cell_idx])

title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = 15, membership_vec = membership_vec,
                                         data_1 = T, data_2 = F,
                                         bool_matrix = T, include_diag = F, verbose = T)

#compute all the degree vectors
k_max <- 50
c_eig <- multiomicCCA::compute_laplacian(rna_frnn$c_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                         colname_vec = paste0("clap_", 1:(k_max-1)))
d_eig <- multiomicCCA::compute_laplacian(rna_frnn$d_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                         colname_vec = paste0("dlap_", 1:(k_max-1)))
e_eig <- multiomicCCA::compute_laplacian(rna_frnn$e_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                         colname_vec = paste0("elap_", 1:(k_max-1)))

################

set.seed(10)
rna_embeddings <- multiomicCCA::plot_embeddings(dcca_res, membership_vec, data_1 = T, data_2 = F, 
                                            add_noise = F, pca = F, only_embedding = T, verbose = T)
mbrain2[["common"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[1]], key = "UMAP", assay = "RNA")
mbrain2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[2]], key = "UMAP", assay = "RNA")
mbrain2[["everything"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[3]], key = "UMAP", assay = "RNA")


########################

p1 <- ncol(mat_1); p2 <- ncol(mat_2)
gene_smoothed <- lapply(1:p1, function(j){
  if(j %% floor(p1/10) == 0) cat('*')
  
  c_res <- compute_smooth_signal(mat_1_denoised[,j], c_eig)
  d_res <- compute_smooth_signal(mat_1_denoised[,j], d_eig)
  e_res <- compute_smooth_signal(mat_1_denoised[,j], e_eig)
  
  list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
       d_variance = d_res$variance, d_r2 = d_res$r_squared,
       e_variance = e_res$variance, e_r2 = e_res$r_squared)
})

# tmp <- sapply(gene_smoothed, function(x){max(x$c_r2, x$d_r2)})
# quantile(tmp)
# tmp <- sapply(gene_smoothed, function(x){x$c_variance})
# round(quantile(tmp),2)
# tmp <- sapply(gene_smoothed, function(x){x$d_variance})
# round(quantile(tmp),2)
# tmp <- sapply(gene_smoothed, function(x){
#   max_val <- max(x$c_variance,x$d_variance)
#   min_val <- min(x$c_variance,x$d_variance)
#   (max_val - min_val)/min_val
# })
# round(quantile(tmp),2)
# length(which(tmp > 0.5))
# sort(colnames(mat_1)[which(tmp > 0.5)])
# idx <- which(colnames(mat_1) %in% gene_vec)
# tmp[idx]

# plot some of the genes with smallest r2
idx <- order(sapply(gene_smoothed, function(x){max(x$c_r2, x$d_r2)}), decreasing = F)[1:5]
for(j in 1:length(idx)){
  c_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], c_eig)
  d_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], d_eig)
  e_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], e_eig)
  
  multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_1)[idx[j]], 
                               prefix = "RNA", e_vec = mat_1_denoised[,idx[j]],
                               c_vec = dcca_decomp$common_mat_1[,idx[j]],
                               d_vec = dcca_decomp$distinct_mat_1[,idx[j]],
                               e_res = e_res, c_res = c_res, d_res = d_res,
                               filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_lowesetr2_gene", j, ".png"))
}

# plot some of the genes with highest r2
idx <- order(sapply(gene_smoothed, function(x){max(x$c_r2, x$d_r2)}), decreasing = T)[1:5]
for(j in 1:length(idx)){
  c_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], c_eig)
  d_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], d_eig)
  e_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], e_eig)
  
  multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_1)[idx[j]], 
                               prefix = "RNA", e_vec = mat_1_denoised[,idx[j]],
                               c_vec = dcca_decomp$common_mat_1[,idx[j]],
                               d_vec = dcca_decomp$distinct_mat_1[,idx[j]],
                               e_res = e_res, c_res = c_res, d_res = d_res,
                               filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_highestr2_gene", j, ".png"))
}

# plot some of the genes with biggest difference
tmp <- sapply(gene_smoothed, function(x){
  max_val <- max(x$c_variance,x$d_variance)
  min_val <- min(x$c_variance,x$d_variance)
  (max_val - min_val)/min_val
})
idx <- order(tmp, decreasing = T)[1:5]
for(j in 1:length(idx)){
  c_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], c_eig)
  d_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], d_eig)
  e_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], e_eig)
  
  multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_1)[idx[j]], 
                               prefix = "RNA", e_vec = mat_1_denoised[,idx[j]],
                               c_vec = dcca_decomp$common_mat_1[,idx[j]],
                               d_vec = dcca_decomp$distinct_mat_1[,idx[j]],
                               e_res = e_res, c_res = c_res, d_res = d_res,
                               filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_biggestdiff_gene", j, ".png"))
}

##########################################
##########################################
##########################################

set.seed(10)
atac_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = 15, membership_vec = membership_vec,
                                         data_1 = F, data_2 = T,
                                         bool_matrix = T, include_diag = F, verbose = T)

#compute all the degree vectors
k_max <- 50
c_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$c_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                         colname_vec = paste0("clap_", 1:(k_max-1)))
d_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$d_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                         colname_vec = paste0("dlap_", 1:(k_max-1)))
e_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$e_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                         colname_vec = paste0("elap_", 1:(k_max-1)))

# atac_smoothed <- lapply(1:p2, function(j){
#   print(j)
#   
#   c_res <- compute_smooth_signal(mat_2_denoised[,j], c_eig2)
#   d_res <- compute_smooth_signal(mat_2_denoised[,j], d_eig2)
#   e_res <- compute_smooth_signal(mat_2_denoised[,j], e_eig2)
#   
#   list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
#        d_variance = d_res$variance, d_r2 = d_res$r_squared,
#        e_variance = e_res$variance, e_r2 = e_res$r_squared)
# })
# 
# save(atac_smoothed, file = "../../../../out/Writeup14b/atac_smoothed.RData")
load("../../../../out/Writeup14b/atac_smoothed.RData")

set.seed(10)
atac_embeddings <- multiomicCCA::plot_embeddings(dcca_res, membership_vec, data_1 = F, data_2 = T, 
                                                add_noise = F, pca = F, only_embedding = T, verbose = T)
mbrain2 <- Seurat::CreateSeuratObject(counts = t(mat_1_denoised))
mbrain2[["label_Savercat"]] <- metadata$label_Savercat[cell_idx]
mbrain2[["common"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[1]], key = "UMAP", assay = "RNA")
mbrain2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[2]], key = "UMAP", assay = "RNA")
mbrain2[["everything"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[3]], key = "UMAP", assay = "RNA")


# plot some of the genes with smallest r2
idx <- order(sapply(atac_smoothed, function(x){max(x$c_r2, x$d_r2)}), decreasing = F)[1:5]
for(j in 1:length(idx)){
  c_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], c_eig)
  d_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], d_eig)
  e_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], e_eig)
  
  multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_2)[idx[j]], 
                               prefix = "ATAC", e_vec = mat_2_denoised[,idx[j]],
                               c_vec = dcca_decomp$common_mat_2[,idx[j]],
                               d_vec = dcca_decomp$distinct_mat_2[,idx[j]],
                               e_res = e_res, c_res = c_res, d_res = d_res,
                               filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_lowesetr2_atac", j, ".png"))
}

# plot some of the genes with highest r2
idx <- order(sapply(atac_smoothed, function(x){max(x$c_r2, x$d_r2)}), decreasing = T)[1:5]
for(j in 1:length(idx)){
  c_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], c_eig)
  d_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], d_eig)
  e_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], e_eig)
  
  multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_2)[idx[j]], 
                               prefix = "ATAC", e_vec = mat_2_denoised[,idx[j]],
                               c_vec = dcca_decomp$common_mat_2[,idx[j]],
                               d_vec = dcca_decomp$distinct_mat_2[,idx[j]],
                               e_res = e_res, c_res = c_res, d_res = d_res,
                               filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_highestr2_atac", j, ".png"))
}

# plot some of the genes with biggest difference
tmp <- sapply(atac_smoothed, function(x){
  max_val <- max(x$c_variance,x$d_variance)
  min_val <- min(x$c_variance,x$d_variance)
  (max_val - min_val)/min_val
})
idx <- order(tmp, decreasing = T)[1:5]
for(j in 1:length(idx)){
  c_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], c_eig)
  d_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], d_eig)
  e_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], e_eig)
  
  multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_2)[idx[j]], 
                               prefix = "ATAC", e_vec = mat_2_denoised[,idx[j]],
                               c_vec = dcca_decomp$common_mat_2[,idx[j]],
                               d_vec = dcca_decomp$distinct_mat_2[,idx[j]],
                               e_res = e_res, c_res = c_res, d_res = d_res,
                               filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_biggestdiff_atac", j, ".png"))
}


########################

png(file = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_laplacian_r2.png"),
    height = 1500, width = 3000, res = 300, units = "px")
par(mfrow = c(1,2))
tmp <- sapply(gene_smoothed, function(x){max(x$c_r2, x$d_r2)})
idx <- which(colnames(mat_1) %in% gene_vec)
col_vec <- rep(rgb(0.5,0.5,0.5,0.2), ncol(mat_1)); col_vec[idx] <- "red"
order_idx <- order(tmp, decreasing = F)
plot(tmp[order_idx], pch = 16, col = col_vec[order_idx], xlab = "Order of genes", ylab = "R2 b/w denoised vector and frNN Laplacian",
     main = "RNA")

idx <- unlist(lapply(1:length(gene_vec), function(i){
  atac_names <- colnames(myobj$X.aggr[[i]]$counts)
  include_bool <- which(myobj$X.aggr[[i]]$peaks.gr$pval.spearman < 0.05)
  atac_names <- atac_names[include_bool]
  atac_idx <- which(colnames(mat_2) %in% atac_names)
  atac_idx
}))
tmp <- sapply(atac_smoothed, function(x){max(x$c_r2, x$d_r2)})
col_vec <- rep(rgb(0.5,0.5,0.5,0.2), ncol(mat_1)); col_vec[idx] <- "red"
order_idx <- order(tmp, decreasing = F)
plot(tmp[order_idx], pch = 16, col = col_vec[order_idx], xlab = "Order of peaks", ylab = "R2 b/w denoised vector and frNN Laplacian",
     main = "ATAC")
graphics.off()

##################



