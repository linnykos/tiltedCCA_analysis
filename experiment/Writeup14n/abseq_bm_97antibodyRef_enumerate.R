rm(list=ls())
load("../../../../out/main/abseq_bm97Ref_preprocessed.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(bm)

Seurat::DefaultAssay(bm) <- "RNA"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:30)
bm <- Seurat::FindClusters(bm, resolution = 0.25)

Seurat::DefaultAssay(bm) <- "AB"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:30, reduction = "apca")
bm <- Seurat::FindClusters(bm, resolution = 0.25)

Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@scale.data[Seurat::VariableFeatures(object = bm),])
Seurat::DefaultAssay(bm) <- "AB"
mat_2 <- Matrix::t(bm[["AB"]]@scale.data)

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

set.seed(10)
multiSVD_obj <- tiltedCCA:::create_multiSVD(mat_1 = mat_1b, mat_2 = mat_2b,
                                            dims_1 = 1:30, dims_2 = 1:30,
                                            center_1 = T, center_2 = T,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = F,
                                            rescale_1 = F, rescale_2 = F,
                                            scale_1 = T, scale_2 = T,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = as.factor(bm$RNA_snn_res.0.25), 
                                           large_clustering_2 = as.factor(bm$AB_snn_res.0.25), 
                                           num_metacells = 5000,
                                           verbose = 1)

#########################

param_df <- data.frame(latent_k =       c(20, 20, 50, 50, 20, 20, 40, 20, 20, 
                                          50, 50, 50, 50,
                                          20, 40, 20,  20,  20,  20, 
                                          20),
                       num_neigh =      c(15, 30, 15, 30, 15, 30, 15, 30, 60,
                                          15, 30, 30, 60,
                                          60, 60, 100, 100, 100, 60,
                                          30),
                       bool_cosine =    c(T,  T,  T,  T,  T,  T,  T,  T,  T,
                                          T,  T,  T,  T,
                                          T,  T,  T,   T,   T,   T,
                                          T,),
                       bool_intersect = c(F,  F,  F,  F,  T,  T,  F,  T,  T,
                                          T,  T,  F,  F,
                                          T,  T,  T,   T,   T,   T,
                                          T),
                       min_deg =        c(15, 30, 15, 30, 15, 30, 15, 15, 30,
                                          15, 30, 15, 30,
                                          15, 15, 100, 15,  30,  60,
                                          30))

for(j in 1:nrow(param_df)){
  tmp_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                      latent_k = param_df$latent_k[j],
                                      num_neigh = param_df$num_neigh[j],
                                      bool_cosine = param_df$bool_cosine[j],
                                      bool_intersect = param_df$bool_intersect[j],
                                      min_deg = param_df$min_deg[j],
                                      verbose = 2)
  
  tmp <- bm; tmp_mat <- tmp_obj$laplacian_list$common_laplacian
  colnames(tmp_mat) <- paste0("tmp_", 1:ncol(tmp_mat))
  set.seed(10); tmp_umap <- Seurat::RunUMAP(tmp_mat)@cell.embeddings
  tmp_umap_full <- matrix(NA, nrow = ncol(tmp), ncol = 2)
  for(i in 1:length(tmp_obj$metacell_obj$metacell_clustering_list)){
    idx <- tmp_obj$metacell_obj$metacell_clustering_list[[i]]
    tmp_umap_full[idx,] <- rep(tmp_umap[i,], each = length(idx))
  }
  set.seed(10)
  tmp_umap_full <- jitter(tmp_umap_full)
  rownames(tmp_umap_full) <- colnames(tmp)
  
  # plot umap of laplacian
  tmp[["common_laplacian"]] <- Seurat::CreateDimReducObject(tmp_umap_full, key = "commonLapUMAP")
  plot1 <- Seurat::DimPlot(tmp, reduction = "common_laplacian",
                           group.by = "celltype", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("k: ", param_df$latent_k[j],
                                           ", nn: ", param_df$num_neigh[j],
                                           ", cos: ", param_df$bool_cosine[j],
                                           "\nint: ", param_df$bool_intersect[j],
                                           ", deg: ", param_df$min_deg[j]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14n/abseq_bm97Ref_laplacian-", j, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  lap_full <- matrix(0, nrow = ncol(tmp), ncol = ncol(tmp_obj$laplacian_list$common_laplacian))
  for(i in 1:length(tmp_obj$metacell_obj$metacell_clustering_list)){
    idx <- tmp_obj$metacell_obj$metacell_clustering_list[[i]]
    lap_full[idx,] <- rep(tmp_obj$laplacian_list$common_laplacian[i,], each = length(idx))
  }
  png(paste0("../../../../out/figures/Writeup14n/abseq_bm97Ref_laplacian_heatmap-", j, ".png"),
      width = 1500, height = 3500, units = "px", res = 300)
  tiltedCCA:::plot_scores_heatmap.list(obj = list(lap_full),
                                       main = paste0("k: ", param_df$latent_k[j],
                                                     ", nn: ", param_df$num_neigh[j],
                                                     ", cos: ", param_df$bool_cosine[j],
                                                     "\nint: ", param_df$bool_intersect[j],
                                                     ", deg: ", param_df$min_deg[j]),
                                       membership_vec = as.factor(tmp$celltype),
                                       log_scale = T, scaling = 2)
  graphics.off()
}



