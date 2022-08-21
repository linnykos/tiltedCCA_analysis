rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)
source("pbmc_224antibody_colorPalette.R")

load("../../../out/main/citeseq_bm25_preprocessed.RData")
load("../../../out/main/citeseq_pbmc224_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(1, ncol(pbmc))
keep_vec[which(pbmc$celltype.l2 == "Doublet")] <- 0
pbmc$keep <- keep_vec
pbmc <- subset(pbmc, keep == 1)

# adjust RNA
# zz <- sapply(pbmc[["SCT"]]@var.features, function(x){which(rownames(pbmc[["SCT"]]) == x)}); zz
Seurat::DefaultAssay(pbmc) <- "SCT"
pbmc[["pca"]] <- NULL
pbmc[["rna.umap"]] <- NULL
pbmc[["SCT"]]@var.features <- pbmc[["SCT"]]@var.features[1:2000]
pbmc[["SCT"]]@scale.data <- pbmc[["SCT"]]@scale.data[pbmc[["SCT"]]@var.features,]
set.seed(10)
pbmc <- Seurat::RunPCA(pbmc, verbose = FALSE)
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:40, assay = 'SCT',
                        reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

# adjust AB
Seurat::DefaultAssay(pbmc) <- "ADT"
pbmc[["apca"]] <- NULL
pbmc[["adt.umap"]] <- NULL
pbmc[["wnn.umap"]] <- NULL
ab_vec <- bm[["ADT"]]@var.features
ab_vec_notin <- ab_vec[which(!ab_vec %in% rownames(pbmc[["ADT"]]))]
ab_vec_notin 
# [1] "CD11a"       "CD127-IL7Ra" "CD197-CCR7"  "CD278-ICOS"  "CD3"        
# [6] "CD38"        "CD4"         "CD56"        "HLA.DR"
sort(rownames(pbmc[["ADT"]]))
ab_vec_other <- c("CD11a/CD18", "CD127", "CD278", "CD3-1", 
                  "CD38-1", "CD38-2", "CD4-2", "CD56-1", 
                  "HLA-DR")
ab_vec_new <- c(ab_vec_other, ab_vec[which(ab_vec %in% rownames(pbmc[["ADT"]]))])
all(ab_vec_new %in% rownames(pbmc[["ADT"]]))
pbmc[["ADT"]]@var.features <- ab_vec_new
pbmc[["ADT"]]@scale.data <- pbmc[["ADT"]]@scale.data[ab_vec_new,]
set.seed(10)
pbmc <- Seurat::RunPCA(pbmc, npcs = 18, reduction.name = 'apca', verbose = F)
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:18, assay = 'AB',
                        reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

##########################################

Seurat::DefaultAssay(pbmc) <- "SCT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:40, reduction = "pca")
set.seed(10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.01)

Seurat::DefaultAssay(pbmc) <- "ADT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:18, reduction = "apca")
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.05)

Seurat::DefaultAssay(pbmc) <- "SCT"
mat_1 <- Matrix::t(pbmc[["SCT"]]@scale.data)
Seurat::DefaultAssay(pbmc) <- "ADT"
mat_2 <- Matrix::t(pbmc[["ADT"]]@scale.data)

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
                                            dims_1 = 1:40, dims_2 = 1:18,
                                            center_1 = T, center_2 = T,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = F,
                                            rescale_1 = F, rescale_2 = F,
                                            scale_1 = T, scale_2 = T,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = as.factor(pbmc$SCT_snn_res.0.01), 
                                           large_clustering_2 = as.factor(pbmc$ADT_snn_res.0.05), 
                                           num_metacells = 5000,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 30,
                                         num_neigh = 30,
                                         bool_cosine = T,
                                         bool_intersect = F,
                                         min_deg = 30,
                                         verbose = 2)

#############

multiSVD_obj <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                      verbose = 1)

multiSVD_obj <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 1)
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = T,
                                                    bool_modality_2_full = T)

save(multiSVD_obj, pbmc,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224-25antibodies_tiltedcca.RData")

#################################

set.seed(10)
pbmc[["common_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                      what = "common",
                                                      aligned_umap_assay = "rna.umap",
                                                      seurat_obj = pbmc,
                                                      seurat_assay = "SCT",
                                                      verbose = 1)
set.seed(10)
pbmc[["distinct1_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                         what = "distinct_1",
                                                         aligned_umap_assay = "rna.umap",
                                                         seurat_obj = pbmc,
                                                         seurat_assay = "SCT",
                                                         verbose = 1)
set.seed(10)
pbmc[["distinct2_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                         what = "distinct_2",
                                                         aligned_umap_assay = "rna.umap",
                                                         seurat_obj = pbmc,
                                                         seurat_assay = "SCT",
                                                         verbose = 1)

save(multiSVD_obj, pbmc,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224-25antibodies_tiltedcca.RData")

