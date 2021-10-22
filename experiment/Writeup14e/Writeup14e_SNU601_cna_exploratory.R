rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/Writeup14e/Writeup14e_SNU_atac_exploratory.RData")

length(which(SNU$nucleosome_signal > 2.75))
ncol(SNU)

SNU[["atac"]]@scale.data <- SNU[["atac"]]@data
SNU[["atac"]]@data <- SNU[["atac"]]@counts
SNU <- subset(SNU, nucleosome_signal < 2.75)

###########################

.cluster_values <- function(mat_rho, mat_theta, seurat_obj){
  rho_avg_list <- lapply(1:max(seurat_obj$clone, na.rm = T), function(k){
    cell_idx <- which(seurat_obj$clone == k)
    vec <- matrixStats::colMeans2(mat_rho[cell_idx,], na.rm = T)
    vec[is.na(vec)] <- 1
    vec
  })
  
  theta_avg_list <- lapply(1:max(seurat_obj$clone, na.rm = T), function(k){
    cell_idx <- which(seurat_obj$clone == k)
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
  
  list(mat_rho = mat_rho, mat_theta = mat_theta)
}

.form_matrix <- function(mat_rho, mat_theta, 
                         rho_const = 1, theta_const = 0.5){
  if(!is.na(rho_const)){
    mat_rho[is.na(mat_rho)] <- rho_const
  } else {
    for(j in 1:ncol(mat_rho)){
      mat_rho[which(is.na(mat_rho[,j])),j] <- mean(mat_rho[,j], na.rm = T)
    }
  }
  
  if(!is.na(theta_const)){
    mat_theta[is.na(mat_theta)] <- theta_const
  } else {
    for(j in 1:ncol(mat_theta)){
      mat_theta[which(is.na(mat_theta[,j])),j] <- mean(mat_theta[,j], na.rm = T)
    }
  }
  
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
  
  mat_2
}

.create_umapobj <- function(mat_2, seurat_obj, rank_k = 7, verbose = T){
  set.seed(10)
  svd_res <- irlba::irlba(mat_2, nv = 50)
  if(verbose) print(svd_res$d)
  dimred <- multiomicCCA:::.mult_mat_vec(svd_res$u[,1:rank_k], svd_res$d[1:rank_k])
  rownames(dimred) <- colnames(seurat_obj)
  colnames(dimred) <- paste0("pca_", 1:ncol(dimred))
  
  set.seed(10)
  umap_res <- Seurat::RunUMAP(dimred)
  umap_res <- umap_res@cell.embeddings
  rownames(umap_res) <- colnames(seurat_obj)
  Seurat::CreateDimReducObject(embedding = umap_res, key = "cnaUMAP_", assay = "atac")
}

###############

genotype_values <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_values_h1_h2_c20.rds")

mat_rho_org <- genotype_values[,grep("^rho*", colnames(genotype_values))]
mat_rho_org <- mat_rho_org[rownames(mat_rho_org) %in% rownames(SNU@meta.data),]
mat_theta_org <- genotype_values[,grep("^theta*", colnames(genotype_values))]
mat_theta_org <- mat_theta_org[rownames(mat_theta_org) %in% rownames(SNU@meta.data),]
all(rownames(mat_rho_org) == rownames(SNU@meta.data))
length(which(is.na(mat_rho_org)))
length(which(is.na(mat_theta_org)))

#################
# version 1: no clustering
mat_2 <- .form_matrix(mat_rho_org, mat_theta_org)
SNU[["cloneumap1"]] <- .create_umapobj(mat_2, SNU, rank_k = 10)

# version 2: clustering
res <- .cluster_values(mat_rho_org, mat_theta_org, SNU)
mat_2 <- .form_matrix(res$mat_rho, res$mat_theta)
SNU[["cloneumap2"]] <- .create_umapobj(mat_2, SNU, rank_k = 10)

##############################

genotype_est <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_est.rds")
canonicalPoints <- data.frame(code = c(1:27),
                              rho = c(0.5,0.5, 1,1,1,1.5,1.5,1.5,1.5,2,2,2,2,2,2.5,2.5,2.5,2.5,2.5,2.5,3,3,3,3,3,3,3),
                              theta = c(0,1,0,0.5,1,0,1/3,2/3,1,0,1/4,2/4,3/4,1,0,1/5,2/5,3/5,4/5,1,0,1/6,2/6,3/6,4/6,5/6,1))
genotype_est <- genotype_est[rownames(genotype_est) %in% rownames(SNU@meta.data),]

mat_rho <- genotype_est
for(i in 1:27){
  mat_rho[genotype_est == i] <- canonicalPoints[i,"rho"]
}
mat_theta <- genotype_est
for(i in 1:27){
  mat_theta[genotype_est == i] <- canonicalPoints[i,"theta"]
}

mat_2 <- .form_matrix(mat_rho, mat_theta)
SNU[["cloneumap3"]] <- .create_umapobj(mat_2, SNU, rank_k = 25)

# version 2: clustering
res <- .cluster_values(mat_rho, mat_theta, SNU)
mat_2 <- .form_matrix(res$mat_rho, res$mat_theta)
SNU[["cloneumap4"]] <- .create_umapobj(mat_2, SNU, rank_k = 25)


###########################

plot1 <- Seurat::DimPlot(SNU, reduction = "cloneumap1",
                         group.by = "clone", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (CNA):\nObserved, no averaging"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

plot2 <- Seurat::DimPlot(SNU, reduction = "cloneumap2",
                         group.by = "clone", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot2 <- plot2 + ggplot2::ggtitle(paste0("SNU (CNA):\nObserved, averaging"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))


plot3 <- Seurat::DimPlot(SNU, reduction = "cloneumap3",
                         group.by = "clone", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot3 <- plot3 + ggplot2::ggtitle(paste0("SNU (CNA):\nEstimated, no averaging"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))


plot4 <- Seurat::DimPlot(SNU, reduction = "cloneumap4",
                         group.by = "clone", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot4 <- plot4 + ggplot2::ggtitle(paste0("SNU (CNA):\nEstimated, averaging"))
plot4 <- plot4 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p <- cowplot::plot_grid(plot1, plot2, plot3, plot4)
cowplot::save_plot(filename = "../../../../out/figures/Writeup14e/Writeup14e_SNU_cna_exploratory_panel.png", p, 
                   ncol = 2, nrow = 2, base_asp = 1.2, device = "png")

