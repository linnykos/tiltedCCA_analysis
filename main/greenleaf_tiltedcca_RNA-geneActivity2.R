rm(list=ls())
load("../../../out/main/10x_greenleaf_preprocessed.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(greenleaf)

Seurat::DefaultAssay(greenleaf) <- "SCT"

var_features <- greenleaf[["SCT"]]@var.features
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
var_features <- unique(c(var_features, mentioned_genes))
var_features <- intersect(var_features, rownames(greenleaf[["SCT"]]))
var_features <- var_features[which(paste0("ATAC-",var_features) %in% rownames(greenleaf[["geneActivity"]]))]

Seurat::DefaultAssay(greenleaf) <- "SCT"
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[var_features,])
Seurat::DefaultAssay(greenleaf) <- "geneActivity"
mat_2 <- Matrix::t(greenleaf[["geneActivity"]]@data[paste0("ATAC-",var_features),])

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
                                            center_1 = T, center_2 = T,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = F,
                                            rescale_1 = F, rescale_2 = F,
                                            scale_1 = T, scale_2 = T,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = NULL, 
                                           large_clustering_2 = NULL, 
                                           num_metacells = 5000,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 20,
                                         num_neigh = 15,
                                         bool_cosine = T,
                                         bool_intersect = F,
                                         min_deg = 15,
                                         verbose = 2)

tmp <- greenleaf; tmp_mat <- multiSVD_obj$laplacian_list$common_laplacian
colnames(tmp_mat) <- paste0("tmp_", 1:ncol(tmp_mat))
set.seed(10); tmp_umap <- Seurat::RunUMAP(tmp_mat)@cell.embeddings
tmp_umap_full <- matrix(NA, nrow = ncol(tmp), ncol = 2)
for(i in 1:length(multiSVD_obj$metacell_obj$metacell_clustering_list)){
  idx <- multiSVD_obj$metacell_obj$metacell_clustering_list[[i]]
  tmp_umap_full[idx,] <- rep(tmp_umap[i,], each = length(idx))
}
set.seed(10)
tmp_umap_full <- jitter(tmp_umap_full)
rownames(tmp_umap_full) <- colnames(tmp)
tmp[["common_laplacian"]] <- Seurat::CreateDimReducObject(tmp_umap_full, key = "commonLapUMAP")
# tmp[["common_laplacian"]] <- Seurat::CreateDimReducObject(tmp_umap, key = "commonLapUMAP")
plot1 <- Seurat::DimPlot(tmp, reduction = "common_laplacian",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+GeneActivity)\nCommon Laplacian"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_umap_RNA-geneActivity2_common-laplacian.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

multiSVD_obj <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                      enforce_boundary = T,
                                      verbose = 1)
multiSVD_obj <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 1)
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = T,
                                                    bool_modality_2_full = T)

save(multiSVD_obj, greenleaf,
     date_of_run, session_info,
     file = "../../../out/main/10x_greenleaf_tcca_RNA-geneActivity2.RData")

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
     file = "../../../out/main/10x_greenleaf_tcca_RNA-geneActivity2.RData")


