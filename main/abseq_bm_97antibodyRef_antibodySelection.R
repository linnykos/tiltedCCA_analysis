rm(list=ls())

library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/abseq_bm97Ref_differential.RData")
load("../../../out/main/abseq_bm97Ref_tcca.RData")
source("bm_97antibodyRef_colorPalette.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(bm) <- "AB"
antibody_names <- rownames(bm)
logpval_vec <- sapply(1:length(antibody_names), function(k){
  if(k %% floor(length(antibody_names)/10) == 0) cat('*')
  antibody <- antibody_names[k]
  
  # cycle through all the celltypes
  celltype_vec <- sapply(1:length(adt_de_list$level_vec), function(i){
    idx <- which(adt_de_list$combn_mat == i, arr.ind = T)[,2]
    vec <-  sapply(idx, function(j){
      idx <- which(rownames(adt_de_list$de_list[[j]]) == antibody)
      if(length(idx) == 0) return(1)
      adt_de_list$de_list[[j]][idx, "p_val"]
    })
    stats::quantile(vec, probs = 0.75)
  })
  
  max(-log10(celltype_vec))
})
logpval_vec <- pmin(logpval_vec, 300)
names(logpval_vec) <- rownames(bm)

Seurat::DefaultAssay(bm) <- "AB"
mat_2 <- Matrix::t(bm[["AB"]]@scale.data)
mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

set.seed(10)
variable_selection_res <- tiltedCCA:::postprocess_distinct_variable_selection(
  input_obj = multiSVD_obj,
  input_mat = mat_2b,
  logpval_vec = logpval_vec,
  cor_threshold = 0.85,
  input_assay = 2,
  min_subsample_cell = NULL,
  seurat_celltype_variable = NULL,
  num_variables = 10,
  seurat_obj = bm,
  verbose = 2
)

# quantile(logpval_vec)
# logpval_vec[variable_selection_res$selected_variables]
# distinct_mat <- multiSVD_obj$distinct_mat_2
# cor_mat <- stats::cor(distinct_mat[,variable_selection_res$selected_variables])

adt_mat2 <- bm[["AB"]]@scale.data
adt_mat2 <- adt_mat2[variable_selection_res$selected_variables[1:5],]
bm[["AB2"]] <- Seurat::CreateAssayObject(counts = adt_mat2)
bm[["AB2"]]@data <- adt_mat2
bm[["AB2"]]@scale.data <- adt_mat2
bm[["AB2"]]@var.features <- rownames(adt_mat2)

Seurat::DefaultAssay(bm) <- "AB2"
bm <- Seurat::RunPCA(bm, reduction.name = 'apca2', 
                     npcs = nrow(adt_mat2), 
                     verbose = F)

set.seed(10)
bm <- Seurat::RunUMAP(bm, reduction = 'apca2', dims = 1:(ncol(bm[["apca2"]]@cell.embeddings)-1), 
                      assay = 'AB2',
                      reduction.name = 'AB2.umap', reduction.key = 'AB2UMAP_')

svd_1 <- multiSVD_obj$svd_1
adt_mat2 <- t(adt_mat2)
consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = NULL, mat_2 = adt_mat2,
                                           dims_1 = NULL, dims_2 = 1:ncol(adt_mat2),
                                           dims_consensus = 1:20,
                                           svd_1 = svd_1, verbose = 1)

bm[["consensusPCA"]] <- Seurat::CreateDimReducObject(consensus_pca$dimred_consensus, 
                                                     assay = "RNA",
                                                     key = "cPC")

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(bm)
bm[["consensusUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                      assay = "RNA",
                                                      key = "cUMAP")

save(variable_selection_res, bm, 
     logpval_vec, 
     consensus_pca,
     date_of_run, session_info,
     file = "../../../out/main/abseq_bm97Ref_varSelect.RData")

variable_selection_res$cor_vec_intial[order(variable_selection_res$logpval_vec, decreasing = T)[1:10]]
variable_selection_res$selected_variables

######################3

plot2 <- Seurat::DimPlot(bm, reduction = "AB2.umap",
                         group.by = "ct", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nADT: Selected antibodies, Thres: ", variable_selection_res$cor_threshold))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_varSelect_adt-umap.png"),
                plot2, device = "png", width = 11, height = 5, units = "in")

plot3 <- Seurat::DimPlot(bm, reduction = "consensusUMAP",
                         group.by = "ct", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nConsensusPCA: Selected antibodies, Thres: ", variable_selection_res$cor_threshold))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_varSelect_consensusPCA-umap.png"),
                plot3, device = "png", width = 11, height = 5, units = "in")
