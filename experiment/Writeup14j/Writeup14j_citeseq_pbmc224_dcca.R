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
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'apca', dims = 1:25, assay = 'ADT',
                        reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

############

n <- ncol(pbmc)
Seurat::DefaultAssay(pbmc) <- "SCT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:40)
set.seed(10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)
tab_vec <- table(pbmc$SCT_snn_res.0.25)
round(tab_vec/n, 3)
rm_idx <- names(tab_vec[tab_vec/n < 0.02])
pbmc$SCT_snn_res.0.25[pbmc$SCT_snn_res.0.25 %in% rm_idx] <- NA

Seurat::DefaultAssay(pbmc) <- "ADT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:50, reduction = "apca")
set.seed(10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)
tab_vec <- table(pbmc$ADT_snn_res.0.25)
round(tab_vec/n, 3)
rm_idx <- names(tab_vec[tab_vec/n < 0.02])
pbmc$ADT_snn_res.0.25[pbmc$ADT_snn_res.0.25 %in% rm_idx] <- NA

############

# make some plots
plot1 <- Seurat::DimPlot(pbmc, reduction = "rna.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = F)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC\nCITE-Seq (RNA, Subset)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_rna_umap_subset.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = F)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC\nCITE-Seq (Protein, Subset)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_adt_umap_subset.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(pbmc, reduction = "rna.umap",
                         group.by = "SCT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = F)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC\nCITE-Seq (RNA, Clustering)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_rna_umap_clustering.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "ADT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = F)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC\nCITE-Seq (Protein, Clustering)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_adt_umap_clustering.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

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

#########################

rank_1 <- 40; rank_2 <- 50; nn <- 30
set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1b, mat_2b, 
                                      dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      num_neigh = nn, 
                                      metacell_clustering_1 = factor(pbmc$SCT_snn_res.0.25),
                                      metacell_clustering_2 = factor(pbmc$ADT_snn_res.0.25),
                                      fix_tilt_perc = F, 
                                      enforce_boundary = F,
                                      verbose = T)
dcca_res$cca_obj
dcca_res$df_percentage
dcca_res$tilt_perc

save(pbmc, dcca_res, 
     rank_1, rank_2, nn, date_of_run, session_info,
     file = "../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca.RData")

dcca_res2 <- multiomicCCA:::fine_tuning(dcca_res, verbose = T)
dcca_res2$tilt_perc

save(pbmc, dcca_res, dcca_res2, 
     rank_1, rank_2, nn, date_of_run, session_info,
     file = "../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca.RData")

