rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(bm)

Seurat::DefaultAssay(bm) <- "RNA"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:30)
bm <- Seurat::FindClusters(bm, resolution = 0.25)
tab_vec <- table(bm$RNA_snn_res.0.25)
round(tab_vec/n, 2)
rm_idx <- names(tab_vec[tab_vec/n < 0.02])
bm$RNA_snn_res.0.25[bm$RNA_snn_res.0.25 %in% rm_idx] <- NA

Seurat::DefaultAssay(bm) <- "ADT"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:18, reduction = "apca")
bm <- Seurat::FindClusters(bm, resolution = 0.25)
tab_vec <- table(bm$ADT_snn_res.0.25)
round(tab_vec/n, 2)
rm_idx <- names(tab_vec[tab_vec/n < 0.02])
bm$ADT_snn_res.0.25[bm$ADT_snn_res.0.25 %in% rm_idx] <- NA

###############

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)

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
                                      metacell_clustering_1 = factor(bm$RNA_snn_res.0.25),
                                      metacell_clustering_2 = factor(bm$ADT_snn_res.0.25),
                                      fix_tilt_perc = F, verbose = T)
dcca_res$cca_obj
dcca_res$df_percentage
dcca_res$tilt_perc

dcca_res2 <- multiomicCCA:::fine_tuning(dcca_res, verbose = T)
dcca_res2$tilt_perc

save(bm, dcca_res, dcca_res2, 
     rank_1, rank_2, nn, date_of_run, session_info,
     file = "../../../../out/Writeup14f/Writeup14f_citeseq_bm_dcca.RData")

###############

# first plot according to celltype
reduction_vec <- c("rna.umap", "adt.umap", "wnn.umap")
group_vec <- c("celltype.l2")
main_vec <- c("(RNA)", "(ATAC)", "(WNN)")
file_vec <- c("rna", "adt", "wnn")

for(i in 1:length(reduction_vec)){
  for(j in 1:length(group_vec)){
    plot1 <- Seurat::DimPlot(bm, reduction = reduction_vec[i],
                             group.by = group_vec[j], label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq)\n", main_vec[i], ": ", group_vec[j]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_", file_vec[i], "_", group_vec[j], ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}

plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap",
                         group.by = "RNA_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq)\nRNA: meta-clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_rna_metacluster.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                         group.by = "ADT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq)\nADT: meta-clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_adt_metacluster.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")




