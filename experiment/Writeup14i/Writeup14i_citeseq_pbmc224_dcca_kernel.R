rm(list=ls())
load("../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_preprocessed.RData")
library(Seurat)
source("kernel_matrix.R")
source("metacell.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(pbmc)

membership_vec <- as.factor(pbmc@meta.data$celltype.l2)
n <- length(membership_vec)
sort(table(membership_vec))
set.seed(10)
idx <- multiomicCCA::construct_celltype_subsample(membership_vec, 
                                                  min_subsample_cell = 5000)
keep_vec <- rep(0, n)
names(keep_vec) <- rownames(pbmc@meta.data)
keep_vec[idx] <- 1
pbmc$keep <- keep_vec
pbmc <- subset(pbmc, keep == 1)
pbmc[["rna.umap"]] <- NULL
pbmc[["adt.umap"]] <- NULL
pbmc[["wnn.umap"]] <- NULL

set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'pca', dims = 1:40, assay = 'SCT',
                        reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'apca', dims = 1:50, assay = 'ADT',
                        reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

############

n <- ncol(pbmc)
Seurat::DefaultAssay(pbmc) <- "SCT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:40)
set.seed(10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)

Seurat::DefaultAssay(pbmc) <- "ADT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:25, reduction = "apca")
set.seed(10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)

###############

mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ADT"]]@scale.data)

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

########################

clustering_1 <- lapply(levels(pbmc$SCT_snn_res.0.25), function(x){
  which(pbmc$SCT_snn_res.0.25 == x)
})
clustering_1 <- clustering_1[sapply(clustering_1, length) != 0]
clustering_2 <- lapply(levels(pbmc$ADT_snn_res.0.25), function(x){
  which(pbmc$ADT_snn_res.0.25 == x)
})
clustering_2 <- clustering_2[sapply(clustering_2, length) != 0]

metacell_clustering_1 <- form_subclusters(mat = mat_1b, 
                                          clustering = clustering_1, 
                                          target_k = 120)
length(which(sapply(metacell_clustering_1, length) > 5))
metacell_clustering_2 <- form_subclusters(mat = mat_2b, 
                                          clustering = clustering_2, 
                                          target_k = 120)
length(which(sapply(metacell_clustering_2, length) > 5))
tmp <- intersect_metacells(metacell_clustering_1 = metacell_clustering_1,
                           metacell_clustering_2 = metacell_clustering_2,
                           n = ncol(pbmc))
metacell_clustering <- tmp$metacell_clustering
clustering_hierarchy_1 <- tmp$clustering_hierarchy_1
clustering_hierarchy_2 <- tmp$clustering_hierarchy_2

#########################

kernel_mat_1 <- form_kernel_matrix(mat = mat_1, K = 40, 
                                   metacell_clustering = metacell_clustering, 
                                   clustering_hierarchy = clustering_hierarchy_1,
                                   symmetrize_func = "average")
kernel_mat_2 <- form_kernel_matrix(mat = mat_2, K = 50, 
                                   metacell_clustering = metacell_clustering, 
                                   clustering_hierarchy = clustering_hierarchy_2,
                                   symmetrize_func = "average")

min_embedding <- compute_min_embedding(kernel_mat_1 = kernel_mat_1, 
                                       kernel_mat_2 = kernel_mat_2,
                                       entrywise_func = "max",
                                       K = 100)
target_embedding <- min_embedding
for(i in 1:ncol(target_embedding)){
  target_embedding[,i] <- target_embedding[,i]/multiomicCCA:::.l2norm(target_embedding[,i])
}

#############################

