rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(multiomicCCA)

load("../../../../out/Writeup14d/Writeup14d_10x_mouseembryo_seurat.RData")
date_of_run <- Sys.time(); session_info <- devtools::session_info()

Seurat::DefaultAssay(mbrain2) <- "SCT"
mat_1 <- Matrix::t(mbrain2[["SCT"]]@data[Seurat::VariableFeatures(object = mbrain2),])
Seurat::DefaultAssay(mbrain2) <- "ATAC"
mat_2 <- Matrix::t(mbrain2[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain2),])

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

##################

rank_1 <- 50; rank_2 <- 50; nn <- 30
set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, 
                                      dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = F,
                                      scale_1 = T, scale_2 = F,
                                      num_neigh = nn, 
                                      metacell_clustering_1 = NA,
                                      metacell_clustering_2 = NA,
                                      fix_tilt_perc = F, verbose = T)
dcca_res$cca_obj
dcca_res$df_percentage
dcca_res$tilt_perc

dcca_res2 <- multiomicCCA:::fine_tuning(dcca_res, verbose = T)
dcca_res2$tilt_perc

save(mbrain2, dcca_res, dcca_res2, 
     rank_1, rank_2, nn, date_of_run, session_info,
     file = "../../../../out/Writeup14f/Writeup14f_10x_mouseembryo_dcca.RData")

#############################

# first plot according to clones
reduction_vec <- c("umap", "umap.atac", "wnn.umap")
group_vec <- c("label_Savercat")
main_vec <- c("(RNA)", "(ATAC)", "(WNN)")
file_vec <- c("rna", "atac", "wnn")

for(i in 1:length(reduction_vec)){
  for(j in 1:length(group_vec)){
    plot1 <- Seurat::DimPlot(mbrain2, reduction = reduction_vec[i],
                             group.by = group_vec[j], label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (10x) ", main_vec[i], ": ", group_vec[j]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_10x_mouseembryo_", file_vec[i], "_", group_vec[j], ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}


