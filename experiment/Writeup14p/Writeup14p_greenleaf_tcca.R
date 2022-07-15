rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/main/10x_greenleaf_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(greenleaf)

Seurat::DefaultAssay(greenleaf) <- "SCT"
var_genes <- Seurat::VariableFeatures(object = greenleaf)

mentioned_genes <- c("ALDH2", "APOE", "AQP4", "ASCL1", "BHLHE22",
                     "BHLHE40", "C16orf89", "CAV2", "DLX2", "DOK5",
                     "DUSP1", "EOMES", "ETV4", "FOS", "FOXJ1", "GLI3",
                     "HAS2", "HES1", "HES4", "HSPA1A", "HSPA1B",
                     "ID3", "IGFBP7", "JUN", "KIF1A", "LIMCH1",
                     "MBP", "MEF2C", "NEUROD1", "NEUROD2", "NEUROD4",
                     "NEUROD6", "NEUROG1", "NEUROG2", "NFIA", "NFIB", "NFIC",
                     "NHLH1", "NR2F1", "PAX6", "RFX4", "RUNX1",
                     "OLIG1", "OLIG2", "SOX2", "SOX3", "SOX6",
                     "SOX9", "SOX10", "SOX21", "SPARCL1", "SNCB", "TBX",
                     "TNC", "TOP2A", "TRB1", "WNT11")
mentioned_genes <- mentioned_genes[mentioned_genes %in% rownames(greenleaf)]
var_genes <- unique(c(var_genes, mentioned_genes))

mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[var_genes,])
Seurat::DefaultAssay(greenleaf) <- "ATAC"
mat_2 <- Matrix::t(greenleaf[["ATAC"]]@data[Seurat::VariableFeatures(object = greenleaf),])

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
                                            dims_1 = 1:50, dims_2 = 2:50,
                                            center_1 = T, center_2 = F,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = T,
                                            rescale_1 = F, rescale_2 = T,
                                            scale_1 = T, scale_2 = F,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = NULL, 
                                           large_clustering_2 = NULL, 
                                           num_metacells = NULL,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 20,
                                         num_neigh = 15,
                                         bool_cosine = T,
                                         bool_intersect = F,
                                         min_deg = 15,
                                         verbose = 2)
##################

multiSVD_obj <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                      enforce_boundary = F,
                                      verbose = 1)
multiSVD_obj <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 1)
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = T,
                                                    bool_modality_2_full = F)

#################################

set.seed(10)
greenleaf[["common_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                           what = "common",
                                                           aligned_umap_assay = "umap",
                                                           seurat_obj = greenleaf,
                                                           seurat_assay = "RNA",
                                                           verbose = 1)
set.seed(10)
greenleaf[["distinct1_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                              what = "distinct_1",
                                                              aligned_umap_assay = "umap",
                                                              seurat_obj = greenleaf,
                                                              seurat_assay = "RNA",
                                                              verbose = 1)
set.seed(10)
greenleaf[["distinct2_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                              what = "distinct_2",
                                                              aligned_umap_assay = "umap",
                                                              seurat_obj = greenleaf,
                                                              seurat_assay = "RNA",
                                                              verbose = 1)

save(multiSVD_obj, greenleaf,
     date_of_run, session_info,
     file = "../../../../out/Writeup14p/Writeup14p_10x_greenleaf_tcca.RData")

