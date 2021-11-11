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

#################################

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

.create_umapobj <- function(mat_2, seurat_obj, rank_k = 40, verbose = T){
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
  
  list(dimred = dimred, umap = umap_res)
}

#################################

genotype_values <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_values_h1_h2_c20.rds")

mat_rho_org <- genotype_values[,grep("^rho*", colnames(genotype_values))]
mat_rho_org <- mat_rho_org[rownames(mat_rho_org) %in% rownames(SNU@meta.data),]
mat_theta_org <- genotype_values[,grep("^theta*", colnames(genotype_values))]
mat_theta_org <- mat_theta_org[rownames(mat_theta_org) %in% rownames(SNU@meta.data),]
all(rownames(mat_rho_org) == rownames(SNU@meta.data))
length(which(is.na(mat_rho_org)))
length(which(is.na(mat_theta_org)))

rank_1 <- 40
mat_2 <- .form_matrix(mat_rho_org, mat_theta_org,
                      rho_const = NA, theta_const = NA)
SNU[["cna"]] <- Seurat::CreateAssayObject(data = t(mat_2))
tmp <- .create_umapobj(mat_2, SNU, rank_k = rank_1)
SNU[["pca"]] <- Seurat::CreateDimReducObject(embedding = tmp$dimred, 
                                                  key = "pca_", 
                                                  assay = "cna")
SNU[["umap.cna"]] <- Seurat::CreateDimReducObject(embedding = tmp$umap, 
                                                  key = "cnaUMAP_", 
                                                  assay = "cna")

#################

set.seed(10)
SNU <- Seurat::FindMultiModalNeighbors(SNU, 
                                       reduction.list = list("pca", "lsi"), 
                                       dims.list = list(1:40, 2:50))
set.seed(10)
SNU <- Seurat::RunUMAP(SNU, 
                       nn.name = "weighted.nn", 
                       reduction.name = "wnn.umap", 
                       reduction.key = "wnnUMAP_")

save(SNU, file = "../../../../out/Writeup14f/Writeup14f_SNU_preprocessed.RData")

######################

# first plot according to clones
reduction_vec <- c("umap.atac", "umap.cna", "wnn.umap")
group_vec <- c("clone")
main_vec <- c("(ATAC)", "(Copy number)", "(WNN)")
file_vec <- c("atac", "cna", "wnn")

for(i in 1:length(reduction_vec)){
  for(j in 1:length(group_vec)){
    plot1 <- Seurat::DimPlot(SNU, reduction = reduction_vec[i],
                             group.by = group_vec[j], label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU ", main_vec[i], ": ", group_vec[j]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14f/Writeup14f_SNU_", file_vec[i], "_", group_vec[j], ".png"),
                    plot1, device = "png", width = 5, height = 5, units = "in")
  }
}

