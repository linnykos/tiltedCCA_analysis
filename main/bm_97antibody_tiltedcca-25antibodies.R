rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/citeseq_bm25_preprocessed.RData")
ab_vec <- bm[["ADT"]]@var.features
load("../../../out/main/abseq_bm97_preprocessed.RData")
source("bm_97antibody_colorPalette.R")
original_var_features <- bm[["RNA"]]@var.features

# adjust RNA
Seurat::DefaultAssay(bm) <- "RNA"
bm <- Seurat::FindVariableFeatures(bm, selection.method = "vst", nfeatures = 2000)
length(bm[["RNA"]]@var.features)
length(intersect(original_var_features, bm[["RNA"]]@var.features))
bm[["pca"]] <- NULL
bm[["rna.umap"]] <- NULL
bm <- Seurat::RunPCA(bm, verbose = F)
set.seed(10)
bm <- Seurat::RunUMAP(bm, reduction = 'pca', dims = 1:30, assay = 'RNA',
                      reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')


# adjust AB
Seurat::DefaultAssay(bm) <- "AB"
bm[["apca"]] <- NULL
bm[["adt.umap"]] <- NULL
bm[["wnn.umap"]] <- NULL
ab_vec2 <- paste0(ab_vec, "-AB")
ab_vec_notin <- ab_vec2[which(!ab_vec2 %in% rownames(bm[["AB"]]))]
ab_vec_notin 
# [1] "CD127-IL7Ra-AB" "CD161-AB"       "CD197-CCR7-AB"  "CD278-ICOS-AB" 
# [5] "CD45RO-AB"      "CD57-AB"        "CD79b-AB"       "CD8a-AB"       
# [9] "HLA.DR-AB"
sort(rownames(bm[["AB"]]))
ab_vec_other <- c("CD127-AB", "CD197-AB", "CD278-AB",
                  "CD8-AB", "CD45RA-AB")
ab_vec_new <- c(ab_vec_other, ab_vec2[which(ab_vec2 %in% rownames(bm[["AB"]]))])
all(ab_vec_new %in% rownames(bm[["AB"]]))
bm[["AB"]]@var.features <- ab_vec_new
bm[["AB"]]@scale.data <- bm[["AB"]]@scale.data[ab_vec_new,]
set.seed(10)
bm <- Seurat::RunPCA(bm, npcs = 18, reduction.name = 'apca', verbose = F)
set.seed(10)
bm <- Seurat::RunUMAP(bm, dims = 1:18, assay = 'AB',
                      reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                         group.by = "ct", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (10x, RNA+ATAC)\nAB UMAP (2000 genes)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97-25antibodies_adt-umap.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

###########

Seurat::DefaultAssay(bm) <- "RNA"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:30)
bm <- Seurat::FindClusters(bm, resolution = 0.25)

Seurat::DefaultAssay(bm) <- "AB"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:18, reduction = "apca")
bm <- Seurat::FindClusters(bm, resolution = 0.25)

plot1 <-Seurat::DimPlot(bm, reduction = "rna.umap",
                        group.by = "RNA_snn_res.0.25", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nRNA clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97-25antibodies_rna-clustering.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "adt.umap",
                        group.by = "AB_snn_res.0.25", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nADT clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97-25antibodies_adt-clustering.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")

#########

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

#################################

set.seed(10)
multiSVD_obj <- tiltedCCA:::create_multiSVD(mat_1 = mat_1b, mat_2 = mat_2b,
                                            dims_1 = 1:30, dims_2 = 1:18,
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
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 20,
                                         num_neigh = 30,
                                         bool_cosine = T,
                                         bool_intersect = F,
                                         min_deg = 30,
                                         verbose = 2)

multiSVD_obj <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                      verbose = 1)
multiSVD_obj <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 1)
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1)

save(multiSVD_obj, bm,
     date_of_run, session_info,
     file = "../../../out/main/abseq_bm97_tcca-25antibodies.RData")

#################################

set.seed(10)
bm[["common_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                    what = "common",
                                                    aligned_umap_assay = "rna.umap",
                                                    seurat_obj = bm,
                                                    seurat_assay = "RNA",
                                                    verbose = 1)
set.seed(10)
bm[["distinct1_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                       what = "distinct_1",
                                                       aligned_umap_assay = "rna.umap",
                                                       seurat_obj = bm,
                                                       seurat_assay = "RNA",
                                                       verbose = 1)
set.seed(10)
bm[["distinct2_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                       what = "distinct_2",
                                                       aligned_umap_assay = "rna.umap",
                                                       seurat_obj = bm,
                                                       seurat_assay = "RNA",
                                                       verbose = 1)

save(multiSVD_obj, bm,
     date_of_run, session_info,
     file = "../../../out/main/abseq_bm97_tcca-25antibodies.RData")

