rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)
library(magrittr)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)
celltype_vec <- bm$celltype.l2

idx <- sample(1:nrow(mat_1), floor(nrow(mat_1)/4))
mat_1 <- mat_1[idx,]
mat_2 <- mat_2[idx,]
celltype_vec <- celltype_vec[idx]

iterations <- 3
cell_fraction <- 1/10
pull_factor <- 0.6

idx_from_list <- list("CD14 Mono", "Memory B", "Naive B")
idx_to_list <- list("CD8 Naive", "CD4 Memory", "CD4 Naive")
for(k in 1:length(idx_from_list)){
  idx_from <- which(celltype_vec %in% idx_from_list[[k]])
  idx_to <- which(celltype_vec %in% idx_to_list[[k]])
  mean_vec_from <- colMeans(mat_2[idx_from,])
  mean_vec_to <- colMeans(mat_2[idx_to,])
  for(i in idx_from){
    mat_2[i,] <- mat_2[i,] - mean_vec_from + mean_vec_to
  }
  
  for(iter in 1:iterations){
    set.seed(iter)
    idx_from <- which(celltype_vec %in% idx_from_list[[k]])
    idx_to <- which(celltype_vec %in% idx_to_list[[k]])
    idx_from_subsample <- sample(idx_from, size = round(length(idx_from)*cell_fraction))
    idx_to_subsample <- sample(idx_to, size = round(length(idx_to)*cell_fraction))
    mean_vec_from <- colMeans(mat_2[idx_from_subsample,])
    mean_vec_to <- colMeans(mat_2[idx_to_subsample,])
    
    for(i in idx_from_subsample){
      mat_2[i,] <- (1-pull_factor) * mat_2[i,] + pull_factor * (mean_vec_from + mean_vec_to)/2
    }
    for(i in idx_to_subsample){
      mat_2[i,] <- (1-pull_factor) * mat_2[i,] + pull_factor * (mean_vec_from + mean_vec_to)/2
    }
    
    all_idx <- c(idx_from_subsample, idx_to_subsample)
    # all_idx <- all_idx[sample(1:length(all_idx), size = floor(length(all_idx)/10))]
    shuf_idx <- sample(all_idx)
    mat_2[all_idx,] <- mat_2[shuf_idx,]
  }
}

idx_from_list <- list(c("Memory B", "Naive B"), c("CD4 Memory", "CD4 Naive"))
idx_to_list <- list("CD14 Mono", "CD8 Naive")
for(k in 1:length(idx_from_list)){
  idx_from <- which(celltype_vec %in% idx_from_list[[k]])
  idx_to <- which(celltype_vec %in% idx_to_list[[k]])
  mean_vec_from <- colMeans(mat_1[idx_from,])
  mean_vec_to <- colMeans(mat_1[idx_to,])
  for(i in idx_from){
    mat_1[i,] <- mat_1[i,] - mean_vec_from + mean_vec_to
  }
  
  for(iter in 1:iterations){
    set.seed(iter)
    idx_from <- which(celltype_vec %in% idx_from_list[[k]])
    idx_to <- which(celltype_vec %in% idx_to_list[[k]])
    idx_from_subsample <-  sample(idx_from, size = round(length(idx_from)*cell_fraction))
    idx_to_subsample <- sample(idx_to, size = round(length(idx_to)*cell_fraction))
    mean_vec_from <- colMeans(mat_1[idx_from_subsample,])
    mean_vec_to <- colMeans(mat_1[idx_to_subsample,])
    
    for(i in idx_from_subsample){
      mat_1[i,] <- (1-pull_factor) * mat_1[i,] + pull_factor * (mean_vec_from + mean_vec_to)/2
    }
    for(i in idx_to_subsample){
      mat_1[i,] <- (1-pull_factor) * mat_1[i,] + pull_factor * (mean_vec_from + mean_vec_to)/2
    }
    
    all_idx <- c(idx_from_subsample, idx_to_subsample)
    # all_idx <- all_idx[sample(1:length(all_idx), size = floor(length(all_idx)/10))]
    shuf_idx <- sample(all_idx)
    mat_1[all_idx,] <- mat_1[shuf_idx,]
  }
}

###############################################

bm2 <- Seurat::CreateSeuratObject(counts = t(mat_1))
bm2[["RNA"]]@data <- t(mat_1)
bm2[["ADT"]] <- Seurat::CreateAssayObject(counts = t(mat_2))
Seurat::DefaultAssay(bm2) <- "ADT"
bm2[["ADT"]]@data <- t(mat_2)
bm2[["ADT"]]@var.features <- colnames(mat_2)
bm2$celltype.l2 <- celltype_vec

keep_vec <- rep(0,ncol(bm2))
keep_vec[which(bm2$celltype.l2 %in% c("CD14 Mono", "Memory B", "Naive B",
                                      "CD8 Naive", "CD4 Memory", "CD4 Naive"))] <- 1
names(keep_vec) <- colnames(bm2)
bm2$keep <- keep_vec
bm2 <- subset(bm2, keep == 1)

Seurat::DefaultAssay(bm2) <- "RNA"
bm2[["RNA"]]@var.features <- colnames(mat_1)
bm2 <- Seurat::ScaleData(bm2)
bm2 <- Seurat::RunPCA(bm2, verbose = F)
Seurat::DefaultAssay(bm2) <- "ADT"
bm2 <- Seurat::ScaleData(bm2)
bm2 <- Seurat::RunPCA(bm2, reduction.name = 'apca', verbose = F)


set.seed(10)
bm2 <- Seurat::RunUMAP(bm2, reduction = 'pca', dims = 1:30, assay = 'RNA',
                       reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

set.seed(10)
bm2 <- Seurat::RunUMAP(bm2, reduction = 'apca', dims = 1:18, assay = 'ADT',
                       reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

set.seed(10)
bm2 <- Seurat::FindMultiModalNeighbors(
  bm2, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)
set.seed(10)
bm2 <- Seurat::RunUMAP(bm2, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                       reduction.key = "wnnUMAP_")

################################

n <- ncol(bm2)
svd_1 <- irlba::irlba(bm2[["RNA"]]@scale.data, 30)
embedding_1 <- multiomicCCA:::.mult_mat_vec(svd_1$v, svd_1$d/svd_1$d[1]*sqrt(n))
svd_2 <- svd(bm2[["ADT"]]@scale.data)
embedding_2 <- multiomicCCA:::.mult_mat_vec(svd_2$v, svd_2$d/svd_2$d[1]*sqrt(n))
embedding_all <- cbind(embedding_1, embedding_2)
rownames(embedding_all) <- colnames(bm2)

embedding_all <- scale(embedding_all, center = T, scale = F)
pca_res <- svd(embedding_all)
consensus_mat <- multiomicCCA:::.mult_mat_vec(pca_res$u[,1:20], pca_res$d[1:20])
rownames(consensus_mat) <- colnames(bm2)

set.seed(10)
consensus_umap <- Seurat::RunUMAP(consensus_mat, 
                                  metric = "euclidean",
                                  reduction.key = "umapConsensusPCA_")
bm2[["consensus.umap"]] <- Seurat::CreateDimReducObject(consensus_umap@cell.embeddings)

################################

anchor_name <- "rna.umap"
other_names <- c("rna.umap", "adt.umap", "wnn.umap", "consensus.umap")

for(umap_name in other_names){
  print(umap_name)
  u_mat1 <- bm[[anchor_name]]@cell.embeddings
  u_mat2 <- bm2[[umap_name]]@cell.embeddings
  u_mat1 <- u_mat1[rownames(u_mat1) %in% rownames(u_mat2),]
  tmp <- svd(t(u_mat1) %*% u_mat2)
  rotation_mat <- tmp$u %*% t(tmp$v)
  tmp <- u_mat2 %*% t(rotation_mat)
  rownames(tmp) <- rownames(bm2@meta.data)
  colnames(tmp) <- colnames(bm2[[umap_name]]@cell.embeddings)
  bm2[[umap_name]]@cell.embeddings <- tmp
}

celltype_table <- c("red", "yellow 1", "yellow 2", "green", "blue 1", "blue")
names(celltype_table) <- sort(unique(bm2$celltype.l2))
# "CD14 Mono", "CD4 Memory", "CD4 Naive", "CD8 Naive", "Memory B", "Naive B"  
color_vec <- scales::hue_pal()(length(unique(bm$celltype.l2)))
names(color_vec) <- sort(unique(bm$celltype.l2))
color_vec <- color_vec[names(color_vec) %in% unique(bm2$celltype.l2)]
names(color_vec) <- celltype_table
celltype_color <- rep(NA, ncol(bm2))
for(i in 1:length(celltype_table)){
  celltype_color[which(bm2$celltype.l2 == names(celltype_table)[i])] <- celltype_table[i]
}
names(celltype_color) <- colnames(bm2)
bm2$celltype_color <- celltype_color

# plot according to clones
reduction_vec <-  c("rna.umap", "adt.umap", "wnn.umap", "consensus.umap")
main_vec <- c("(Modality 1)", "(Modality 2)", "(WNN)", "(Consensus PCA)")
file_vec <- c("modality1", "modality2", "wnn", "consensuspca")

for(i in 1:length(reduction_vec)){
  plot1 <- Seurat::DimPlot(bm2, reduction = reduction_vec[i],
                           group.by = "celltype_color", label = TRUE,
                           repel = TRUE, label.size = 2.5,
                           cols = color_vec)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real:\n", main_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_pseudoreal_", file_vec[i], ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# ######################################
# ######################################
# dims_1 <- 1:30; dims_2 <- 1:18; nn <- 30
# 
# rank_1 <- max(dims_1); rank_2 <- max(dims_2)
# n <- nrow(mat_1)
# 
# svd_1 <- multiomicCCA:::.svd_truncated(mat_1, K = rank_1, symmetric = F, rescale = F, 
#                                        mean_vec = F, sd_vec = F, K_full_rank = F)
# svd_2 <- multiomicCCA:::.svd_truncated(mat_2, K = rank_2, symmetric = F, rescale = F, 
#                                        mean_vec = F, sd_vec = F, K_full_rank = F)
# 
# svd_1 <- multiomicCCA:::.check_svd(svd_1, dims = dims_1)
# svd_2 <- multiomicCCA:::.check_svd(svd_2, dims = dims_2)
# 
# dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
# dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)
# 
# set.seed(10)
# jive_umap <- Seurat::RunUMAP(jive_res$embedding, 
#                              metric = "cosine",
#                              reduction.key = "umapJive1_")
# rownames(jive_umap@cell.embeddings) <- rownames(bm@meta.data)
# bm[["jive_umap"]] <- Seurat::CreateDimReducObject(jive_umap@cell.embeddings)
# 

  