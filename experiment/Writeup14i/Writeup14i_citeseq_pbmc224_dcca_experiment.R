rm(list=ls())
load("../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_preprocessed.RData")
library(Seurat)

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
# tab_vec <- table(pbmc$SCT_snn_res.0.25)
# round(tab_vec/n, 3)
# rm_idx <- names(tab_vec[tab_vec/n < 0.02])
# pbmc$SCT_snn_res.0.25[pbmc$SCT_snn_res.0.25 %in% rm_idx] <- NA

Seurat::DefaultAssay(pbmc) <- "ADT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:25, reduction = "apca")
set.seed(10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)
# tab_vec <- table(pbmc$ADT_snn_res.0.25)
# round(tab_vec/n, 3)
# rm_idx <- names(tab_vec[tab_vec/n < 0.02])
# pbmc$ADT_snn_res.0.25[pbmc$ADT_snn_res.0.25 %in% rm_idx] <- NA

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

# do some checks
for(i in 1:length(clustering_hierarchy_1)){
  metacells <- clustering_hierarchy_1[[i]]
  cell_idx <- unlist(metacell_clustering[metacells])
  stopifnot(length(unique(pbmc$SCT_snn_res.0.25[cell_idx])) == 1)
}
for(i in 1:length(clustering_hierarchy_2)){
  metacells <- clustering_hierarchy_2[[i]]
  cell_idx <- unlist(metacell_clustering[metacells])
  stopifnot(length(unique(pbmc$ADT_snn_res.0.25[cell_idx])) == 1)
}

quantile(sapply(metacell_clustering, length), probs = seq(0,1,length.out=11))
length(metacell_clustering)
length(which(sapply(metacell_clustering, length) >= 5))
###############################

shuff_vec <- sample(1:length(metacell_clustering_1))
metacell_vec_1 <- rep(NA, ncol(pbmc))
for(i in 1:length(metacell_clustering_1)){
  metacell_vec_1[metacell_clustering_1[[i]]] <- shuff_vec[i]
}
shuff_vec <- sample(1:length(metacell_clustering_2))
metacell_vec_2 <- rep(NA, ncol(pbmc))
for(i in 1:length(metacell_clustering_2)){
  metacell_vec_2[metacell_clustering_2[[i]]] <- shuff_vec[i]
}
pbmc$metacell_vec_1 <- metacell_vec_1
pbmc$metacell_vec_2 <- metacell_vec_2

plot1 <- Seurat::DimPlot(pbmc, reduction = "rna.umap",
                         group.by = "metacell_vec_1", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nRNA meta-clusters"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14i/Writeup14i_citeseq_pbmc_rna_metacluster.png"),
                plot1, device = "png", width = 8, height = 5, units = "in")
plot1 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "metacell_vec_2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nADT meta-clusters"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14i/Writeup14i_citeseq_pbmc_adt_metacluster.png"),
                plot1, device = "png", width = 8, height = 5, units = "in")


plot1 <- Seurat::DimPlot(pbmc, reduction = "rna.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nRNA cell-types"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14i/Writeup14i_citeseq_pbmc_rna_celltype.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")
plot1 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nADT cell-types"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14i/Writeup14i_citeseq_pbmc_adt_celltype.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")


plot1 <- Seurat::DimPlot(pbmc, reduction = "rna.umap",
                         group.by = "SCT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nRNA large clusters"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14i/Writeup14i_citeseq_pbmc_rna_cluster.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")
plot1 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "ADT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nADT large clusters"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14i/Writeup14i_citeseq_pbmc_adt_cluster.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


#####################

color_df <- data.frame(celltype = sort(unique(pbmc$celltype.l2)),
                       color = scales::hue_pal()(length(unique(pbmc$celltype.l2))))
metacell_df <- as.data.frame(t(sapply(metacell_clustering, function(vec){
  celltype_tab <- table(pbmc$celltype.l2[vec])
  celltype_name <- names(celltype_tab)[which.max(celltype_tab)]
  col <- color_df[which(color_df$celltype == celltype_name), "color"]
  
  c(celltype = celltype_name, color = col)
})))

rna_avg_umap <- t(sapply(metacell_clustering, function(vec){
  matrixStats::colMeans2(pbmc[["rna.umap"]]@cell.embeddings[vec,,drop = F])
}))
png("../../../../out/figures/Writeup14i/Writeup14i_citeseq_pbmc_rna_avg_metacell.png",
    height = 2500, width = 2500, res = 300, units = "px")
plot(rna_avg_umap[,1], rna_avg_umap[,2], col = metacell_col,
     asp = T, pch = 16, 
     xlab = "rnaUMAP_1", ylab = "rnaUMAP_2",
     main = "Human PBMC (Cite-seq):\nRNA metacell averages")
graphics.off()

adt_avg_umap <- t(sapply(metacell_clustering, function(vec){
  matrixStats::colMeans2(pbmc[["adt.umap"]]@cell.embeddings[vec,,drop = F])
}))
png("../../../../out/figures/Writeup14i/Writeup14i_citeseq_pbmc_adt_avg_metacell.png",
    height = 2500, width = 2500, res = 300, units = "px")
plot(adt_avg_umap[,1], adt_avg_umap[,2], col = metacell_col,
     asp = T, pch = 16, 
     xlab = "adtUMAP_1", ylab = "adtUMAP_2",
     main = "Human PBMC (Cite-seq):\nADT metacell averages")
graphics.off()

#########################3

kernel_mat_1 <- form_kernel_matrix(mat = mat_1, K = 40, 
                                   metacell_clustering = metacell_clustering, 
                                   clustering_hierarchy = clustering_hierarchy_1,
                                   symmetrize_func = "average")
kernel_mat_2 <- form_kernel_matrix(mat = mat_2, K = 50, 
                                   metacell_clustering = metacell_clustering, 
                                   clustering_hierarchy = clustering_hierarchy_2,
                                   symmetrize_func = "average")

# cd4_idx <- grep("CD4", metacell_df$celltype)
# cd8_idx <- grep("CD8", metacell_df$celltype)
# tmp1 <- kernel_mat_1; diag(tmp1) <- NA
# zz1 <- matrix(NA, 2, 2)
# zz1[1,1] <- mean(tmp1[cd4_idx, cd4_idx], na.rm = T)
# zz1[1,2] <- mean(tmp1[cd4_idx, cd8_idx], na.rm = T); zz1[2,1] <- zz1[1,2]
# zz1[2,2] <- mean(tmp1[cd8_idx, cd8_idx], na.rm = T)
# tmp2 <- kernel_mat_2; diag(tmp2) <- NA
# zz2 <- matrix(NA, 2, 2)
# zz2[1,1] <- mean(tmp2[cd4_idx, cd4_idx], na.rm = T)
# zz2[1,2] <- mean(tmp2[cd4_idx, cd8_idx], na.rm = T); zz2[2,1] <- zz2[1,2]
# zz2[2,2] <- mean(tmp2[cd8_idx, cd8_idx], na.rm = T)
# zz1; zz2
idx <- which(metacell_df$celltype == "CD14 Mono")
round(kernel_mat_1[idx,idx],2)
round(kernel_mat_2[idx,idx],2)

#######################3

min_embedding <- compute_min_embedding(kernel_mat_1 = kernel_mat_1, 
                                       kernel_mat_2 = kernel_mat_2,
                                       entrywise_func = "max",
                                       K = 100)

set.seed(10)
min_umap <- Seurat::RunUMAP(min_embedding)@cell.embeddings
png("../../../../out/figures/Writeup14i/Writeup14i_citeseq_pbmc_min_umap.png",
    height = 2500, width = 2500, res = 300, units = "px")
plot(min_umap[,1], min_umap[,2], col = metacell_col,
     asp = T, pch = 16, 
     xlab = "minUMAP_1", ylab = "minUMAP_2",
     main = "Human PBMC (Cite-seq):\nMin-embedding metacell averages")
graphics.off()

# celltype_group_list <- list(
#   Mono = c("ASDC", "cDC1", "cDC2", "CD14 Mono", "CD16 Mono", "Doublet", "pDC",  "Platelet"),
#   B = c("B intermediate", "B memory", "B naive", "Plasmablast"),
#   CD4_CD8 = c("Eryth", "CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM",
#     "CD4 TEM", "CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM", "dnT",
#     "gdT", "HSPC", "ILC", "MAIT", "NK", "NK Proliferating", "NK_CD56bright", "Treg")
# )

max_embedding <- compute_min_embedding(kernel_mat_1 = kernel_mat_1, 
                                       kernel_mat_2 = kernel_mat_2,
                                       entrywise_func = "min",
                                       K = 100)

set.seed(10)
max_umap <- Seurat::RunUMAP(max_embedding)@cell.embeddings
png("../../../../out/figures/Writeup14i/Writeup14i_citeseq_pbmc_max_umap.png",
    height = 2500, width = 2500, res = 300, units = "px")
plot(max_umap[,1], max_umap[,2], col = metacell_col,
     asp = T, pch = 16, 
     xlab = "maxUMAP_1", ylab = "maxUMAP_2",
     main = "Human PBMC (Cite-seq):\nMax-embedding metacell averages")
graphics.off()

