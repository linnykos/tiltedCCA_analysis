rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")

date_of_run <- Sys.time(); session_info <- sessionInfo()

mat_1 <- t(mbrain[["SCT"]]@scale.data)
Seurat::DefaultAssay(mbrain) <- "ATAC"
mat_2 <- Matrix::t(mbrain[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain),])

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
metadata <- mbrain@meta.data

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

##########

set.seed(10)
obj <- multiomicCCA:::.prepare_embeddings_singleton(dcca_res$common_score, dcca_res$distinct_score_1,
                                             dcca_res$svd_1, add_noise = F)
set.seed(10)
seurat_umap <- Seurat::RunUMAP(obj$common, metric = "euclidean")

mbrain2 <- Seurat::CreateSeuratObject(counts = t(mat_1))
mbrain2[["label_Savercat"]] <- metadata$label_Savercat[cell_idx]
mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = seurat_umap@cell.embeddings, 
                                                  key = "UMAP", assay = "RNA")
plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", label = TRUE,
                         group.by = "label_Savercat",
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_umap_original_test.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##################

# now try recreating it via the graph
set.seed(10)
rann_res <- RANN::nn2(obj$common, k = 30)
rann_res$nn.idx <- rann_res$nn.idx[,-1]
rann_res$nn.dist <- rann_res$nn.dist[,-1]

rann_res$nn.idx <- lapply(1:nrow(rann_res$nn.idx), function(i){
  rann_res$nn.idx[i,]
})
rann_res$nn.dist <- lapply(1:nrow(rann_res$nn.dist), function(i){
  rann_res$nn.dist[i,]
})

j_vec <- unlist(rann_res$nn.idx)
i_vec <- unlist(lapply(1:length(rann_res$nn.idx), function(i){
  rep(i, length(rann_res$nn.idx[[i]]))
}))
x_vec <- unlist(rann_res$nn.dist)
dist_mat <- Matrix::sparseMatrix(i = i_vec, j = j_vec, x = x_vec)

rownames(dist_mat) <- rownames(mat_1)
colnames(dist_mat) <- rownames(mat_1)

embedding <- custom_umap(dist_mat)
mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = embedding, 
                                                  key = "UMAP", assay = "RNA")
plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", label = TRUE,
                         group.by = "label_Savercat",
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_umap_original_test4.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


embedding <- custom_umap(dist_mat)
mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = embedding, 
                                                  key = "UMAP", assay = "RNA")
plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", label = TRUE,
                         group.by = "label_Savercat",
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_umap_original_test4b.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")
