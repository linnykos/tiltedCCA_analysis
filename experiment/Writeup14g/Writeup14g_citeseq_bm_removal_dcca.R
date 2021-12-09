rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

bm2 <- bm
Seurat::DefaultAssay(bm2) <- "RNA"
bm2[["ADT"]] <- NULL
bm2[["spca"]] <- NULL; bm2[["adt.umap"]] <- NULL; bm2[["wnn.umap"]] <- NULL

mat <- bm[["ADT"]]@counts
protein_idx <- which(!rownames(mat) %in% c("CD8a", "CD4"))
mat <- mat[protein_idx,]
bm2[["ADT"]] <- Seurat::CreateAssayObject(counts = mat)
Seurat::DefaultAssay(bm2) <- "ADT"
bm2 <- Seurat::NormalizeData(bm2, normalization.method = 'CLR', margin = 2)
bm2 <- Seurat::ScaleData(bm2) 
bm2[["ADT"]]@var.features <- rownames(bm2)
bm2 <- Seurat::RunPCA(bm2, reduction.name = 'apca')

set.seed(10)
bm2 <- Seurat::RunUMAP(bm2, reduction = 'apca', dims = 1:18, assay = 'ADT',
                      reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

u_mat1 <- bm2[["rna.umap"]]@cell.embeddings
u_mat2 <- bm2[["adt.umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(bm2@meta.data)
colnames(tmp) <- colnames(bm2[["adt.umap"]]@cell.embeddings)
bm2[["adt.umap"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(bm2, reduction = "adt.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq)\nADT: Removed CD4 and CD8a"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_adt_removal_celltype.l2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

###########################

Seurat::DefaultAssay(bm2) <- "RNA"
set.seed(10)
bm2 <- Seurat::FindNeighbors(bm2, dims = 1:30)
bm2 <- Seurat::FindClusters(bm2, resolution = 0.25)
tab_vec <- table(bm2$RNA_snn_res.0.25)
round(tab_vec/n, 2)
rm_idx <- names(tab_vec[tab_vec/n < 0.02])
bm2$RNA_snn_res.0.25[bm2$RNA_snn_res.0.25 %in% rm_idx] <- NA

Seurat::DefaultAssay(bm2) <- "ADT"
set.seed(10)
bm2 <- Seurat::FindNeighbors(bm2, dims = 1:18, reduction = "apca")
bm2 <- Seurat::FindClusters(bm2, resolution = 0.25)
tab_vec <- table(bm2$ADT_snn_res.0.25)
round(tab_vec/n, 2)
rm_idx <- names(tab_vec[tab_vec/n < 0.02])
bm2$ADT_snn_res.0.25[bm2$ADT_snn_res.0.25 %in% rm_idx] <- NA

plot1 <- Seurat::DimPlot(bm2, reduction = "rna.umap",
                         group.by = "RNA_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq) RNA:\nMeta-clustering (Removed CD4 and CD8a)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_rna_removal_metacluster.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(bm2, reduction = "adt.umap",
                         group.by = "ADT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq)ADT:\nMeta-clustering (Removed CD4 and CD8a)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_adt_removal_metacluster.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

########################

mat_1 <- t(bm2[["RNA"]]@scale.data)
mat_2 <- t(bm2[["ADT"]]@scale.data)

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

#########################

rank_1 <- 30; rank_2 <- 18; nn <- 30
set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1b, mat_2b, 
                                      dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      num_neigh = nn, 
                                      metacell_clustering_1 = factor(bm2$RNA_snn_res.0.25),
                                      metacell_clustering_2 = factor(bm2$ADT_snn_res.0.25),
                                      fix_tilt_perc = F, verbose = T)
dcca_res$cca_obj
dcca_res$df_percentage
dcca_res$tilt_perc

dcca_res2 <- multiomicCCA:::fine_tuning(dcca_res, verbose = T)
dcca_res2$tilt_perc

save(bm2, dcca_res, dcca_res2, 
     rank_1, rank_2, nn, date_of_run, session_info,
     file = "../../../../out/Writeup14g/Writeup14g_citeseq_bm_removal_dcca.RData")



