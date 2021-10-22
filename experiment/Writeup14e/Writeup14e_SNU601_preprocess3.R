rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/Writeup14e/Writeup14e_SNU_atac_exploratory.RData")

length(which(SNU$nucleosome_signal > 2.75))
ncol(SNU)

SNU[["atac"]]@scale.data <- SNU[["atac"]]@data
SNU[["atac"]]@data <- SNU[["atac"]]@counts
SNU <- subset(SNU, nucleosome_signal < 2.75)

set.seed(10)
SNU <- Seurat::RunUMAP(SNU, reduction = 'lsi', dims = 1:50, 
                       reduction.name = "umap.atac", 
                       reduction.key = "atacUMAP_")

################

genotype_values <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_values_h1_h2_c20.rds")

mat_rho <- genotype_values[,grep("^rho*", colnames(genotype_values))]
mat_rho <- mat_rho[rownames(mat_rho) %in% rownames(SNU@meta.data),]
mat_theta <- genotype_values[,grep("^theta*", colnames(genotype_values))]
mat_theta <- mat_theta[rownames(mat_theta) %in% rownames(SNU@meta.data),]
all(rownames(mat_rho) == rownames(SNU@meta.data))

# resolve NA's
rho_avg_list <- lapply(1:max(SNU$clone, na.rm = T), function(k){
  cell_idx <- which(SNU$clone == k)
  vec <- matrixStats::colMeans2(mat_rho[cell_idx,], na.rm = T)
  vec[is.na(vec)] <- 1
  vec
})
theta_avg_list <- lapply(1:max(SNU$clone, na.rm = T), function(k){
  cell_idx <- which(SNU$clone == k)
  vec <- matrixStats::colMeans2(mat_theta[cell_idx,],na.rm = T)
  vec[is.na(vec)] <- 0.5
  vec
})

for(k in 1:max(SNU$clone, na.rm = T)){
  cell_idx <- which(SNU$clone == k)
  for(j in 1:ncol(mat_rho)){
    mat_rho[cell_idx[which(is.na(mat_rho[cell_idx,j]))], j] <- rho_avg_list[[k]][j]
  }
}

for(k in 1:max(SNU$clone, na.rm = T)){
  cell_idx <- which(SNU$clone == k)
  for(j in 1:ncol(mat_theta)){
    mat_theta[cell_idx[which(is.na(mat_theta[cell_idx,j]))], j] <- theta_avg_list[[k]][j]
  }
}

mat_rho[is.na(mat_rho)] <- 1
mat_theta[is.na(mat_theta)] <- 0.5

mat_rho <- mat_rho - 1
mat_rho <- mat_rho/svd(mat_rho)$d[1]*100
mat_theta <- mat_theta - 0.5
mat_theta <- mat_theta/svd(mat_theta)$d[1]*100
mat_2 <- cbind(mat_rho, mat_theta)
colnames(mat_2)[1:ncol(mat_rho)] <- paste0("rho-",colnames(mat_2)[1:ncol(mat_rho)])
colnames(mat_2)[ncol(mat_rho)+(1:ncol(mat_theta))] <- paste0("theta-",colnames(mat_2)[ncol(mat_rho)+(1:ncol(mat_theta))])

sd_vec <- matrixStats::colSds(mat_2)
if(any(sd_vec < 1e-6)){
  idx <- which(sd_vec < 1e-6)
  mat_2 <- mat_2[,-idx]
}

SNU[["CNA"]] <- Seurat::CreateAssayObject(counts = t(mat_2))
set.seed(10)
svd_res <- irlba::irlba(mat_2, nv = 50)
cna_k <- 7
dimred <- multiomicCCA:::.mult_mat_vec(svd_res$u[,1:cna_k], svd_res$d[1:cna_k])
rownames(dimred) <- colnames(SNU)
colnames(dimred) <- paste0("pca_", 1:ncol(dimred))
SNU[["pca"]] <- Seurat::CreateDimReducObject(embedding = dimred, 
                                             key = "pca_", assay = "CNA")
set.seed(10)
umap_res <- Seurat::RunUMAP(dimred)
umap_res <- umap_res@cell.embeddings
rownames(umap_res) <- colnames(SNU)
SNU[["umap.cna"]] <- Seurat::CreateDimReducObject(embedding = umap_res, 
                                                  key = "cnaUMAP_", assay = "CNA")

set.seed(10)
SNU <- Seurat::FindMultiModalNeighbors(SNU, 
                                       reduction.list = list("pca", "lsi"), 
                                       dims.list = list(1:7, 1:50))
set.seed(10)
SNU <- Seurat::RunUMAP(SNU, 
                       nn.name = "weighted.nn", 
                       reduction.name = "wnn.umap", 
                       reduction.key = "wnnUMAP_")

description <- c("CNA based on observed theta/rhos with 7 latent dimensions. ATAC formed by regressing nCount_atac out of ATAC and then applying Signac::RunSVD. (Notably, not using TF-IDF)")
save(SNU, description,
     file = "../../../../out/Writeup14e/Writeup14e_SNU_preprocess3.RData")

#################3

group_vec <- c("clone")
for(group in group_vec){
  plot1 <- Seurat::DimPlot(SNU, reduction = "umap.atac",
                           group.by = group, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (ATAC)\n",group))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_atac_preprocess3_", group, ".png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
}

group_vec <- c("nCount_atac", "nFeature_atac", "total_frags",
               "frac_peak", "frac_promoter", "frac_tss",
               "frac_enhancer", "nucleosome_signal",
               "nucleosome_percentile", "TSS.enrichment", "TSS.percentile")
for(group in group_vec){
  plot1 <- Seurat::FeaturePlot(SNU, reduction = "umap.atac",
                               features = group)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (ATAC)\n",group))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_atac_preprocess3_", group, ".png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
}

###############

p_list <- lapply(sort(unique(SNU@meta.data[,"clone"])), function(clone){
  p0 <- Seurat::DimPlot(SNU, reduction = "umap.atac",
                        cells.highlight = rownames(SNU@meta.data)[which(SNU[["clone"]] == clone)])
  p0 <- p0 + ggplot2::theme(legend.position="none") + ggplot2::ggtitle(paste0("Clone ", clone))
  p0
})

p <- cowplot::plot_grid(p_list[[1]], p_list[[2]], p_list[[3]], 
                        p_list[[4]], p_list[[5]], p_list[[6]])
cowplot::save_plot(filename = "../../../../out/figures/Writeup14e/Writeup14e_SNU_atac_preprocess3_clone_separate.png", p, 
                   ncol = 3, nrow = 2, base_asp = 1, device = "png")

##################

plot1 <- Seurat::DimPlot(SNU, reduction = "umap.cna",
                         group.by = "clone", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (CNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_cna_preprocess3.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(SNU, reduction = "wnn.umap",
                         group.by = "clone", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (CNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_wnn_preprocess3.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <-  Seurat::VlnPlot(SNU, features = "CNA.weight", 
                          group.by = "clone", sort = TRUE, 
                          pt.size = 0.1) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (CNA): WNN weights\nAverage: ", round(SNU$CNA.weight)))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_wnn_preprocess3_weights.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

###############

p_list <- lapply(sort(unique(SNU@meta.data[,"clone"])), function(clone){
  p0 <- Seurat::DimPlot(SNU, reduction = "umap.cna",
                        cells.highlight = rownames(SNU@meta.data)[which(SNU[["clone"]] == clone)])
  p0 <- p0 + ggplot2::theme(legend.position="none") + ggplot2::ggtitle(paste0("Clone ", clone))
  p0
})

p <- cowplot::plot_grid(p_list[[1]], p_list[[2]], p_list[[3]], 
                        p_list[[4]], p_list[[5]], p_list[[6]])
cowplot::save_plot(filename = "../../../../out/figures/Writeup14e/Writeup14e_SNU_cna_preprocess3_clone_separate.png", p, 
                   ncol = 3, nrow = 2, base_asp = 1, device = "png")
