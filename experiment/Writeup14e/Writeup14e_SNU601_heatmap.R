rm(list=ls())
load("../../../../out/Writeup14e/Writeup14e_SNU_atac_exploratory.RData")

length(which(SNU$nucleosome_signal > 2.75))
ncol(SNU)

SNU[["atac"]]@scale.data <- SNU[["atac"]]@data
SNU[["atac"]]@data <- SNU[["atac"]]@counts
SNU <- subset(SNU, nucleosome_signal < 2.75)

genotype_values <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_values_h1_h2_c20.rds")

#############

col_vec <- c(rgb(108,152,208, maxColorValue = 255), #-m
             rgb(46,43,141, maxColorValue = 255), #M-
             rgb(169,205,226, maxColorValue = 255), #mm
             rgb(177,172,173, maxColorValue = 255), #Mm
             rgb(117,198,227, maxColorValue = 255), #MM
             rgb(254,242,186, maxColorValue = 255), #mmm
             rgb(254,171,20, maxColorValue = 255), #Mmm
             rgb(214,102,32, maxColorValue = 255), #MMm
             rgb(177,60,22, maxColorValue = 255), #MMM
             rgb(254,233,231, maxColorValue = 255), #mmmm
             rgb(255,205,199, maxColorValue = 255), #Mmmm
             rgb(236,92,155, maxColorValue = 255), #MMmm
             rgb(171,25,131, maxColorValue = 255), #MMMm
             rgb(69,24,102, maxColorValue = 255), #MMMM
             rgb(222,247,187, maxColorValue = 255), #mmmmm
             rgb(151,231,212, maxColorValue = 255), #Mmmmm
             rgb(87,188,203, maxColorValue = 255), #MMmmm
             rgb(63,169,71, maxColorValue = 255), #MMMmm
             rgb(16,105,44, maxColorValue = 255), #MMMMm
             rgb(18,20,19, maxColorValue = 255), #MMMMMM
             rgb(116,37,26, maxColorValue = 255)) #+
             
mat_rho <- genotype_values[,grep("^rho*", colnames(genotype_values))]
mat_rho <- mat_rho[rownames(mat_rho) %in% rownames(SNU@meta.data),]
mat_theta <- genotype_values[,grep("^theta*", colnames(genotype_values))]
mat_theta <- mat_theta[rownames(mat_theta) %in% rownames(SNU@meta.data),]
mat_rho[is.na(mat_rho)] <- 1
mat_theta[is.na(mat_theta)] <- 0.5

#############################

color_mat <- matrix(NA, nrow = nrow(mat_theta), ncol = ncol(mat_theta))
for(i in 1:nrow(mat_theta)){
  for(j in 1:ncol(mat_theta)){
    rho <- round(mat_rho[i,j]*2)/2
    if(rho == 0){
      color_mat[i,j] <- 0
    } else if(rho == 0.5){
      theta <- round(mat_rho[i,j])
      color_mat[i,j] <- ifelse(theta == 0, 1, 2) 
    
    } else if(abs(rho - 1) <= 1e-3){
      theta <- round(mat_theta[i,j]*2)/2
      if(abs(theta) <= 1e-3){
        color_mat[i,j] <- 3
      } else if(abs(theta - 0.5) <= 1e-3){
        color_mat[i,j] <- 4
      } else {
        color_mat[i,j] <- 5
      }
      
    } else if(abs(rho - 1.5) <= 1e-3){
      theta <- round(mat_theta[i,j]*3)/3
      if((abs(theta) <= 1e-3)){
        color_mat[i,j] <- 6
      } else if(abs(theta - 1/3) <= 1e-3){
        color_mat[i,j] <- 7
      } else if(abs(theta - 2/3) <= 1e-3){
        color_mat[i,j] <- 8
      } else {
        color_mat[i,j] <- 9
      }
      
    } else if(abs(rho - 2) <= 1e-3){
      theta <- round(mat_theta[i,j]*4)/4
      if(abs(theta) <= 1e-3){
        color_mat[i,j] <- 10
      } else if(abs(theta - 1/4) <= 1e-3){
        color_mat[i,j] <- 11
      } else if(abs(theta - 2/4) <= 1e-3){
        color_mat[i,j] <- 12
      } else if(abs(theta - 3/4) <= 1e-3){
        color_mat[i,j] <- 13
      } else {
        color_mat[i,j] <- 14
      }
      
    } else if(abs(rho - 2.5) <= 1e-3){
      theta <- round(mat_theta[i,j]*5)/5
      if(abs(theta) <= 1e-3){
        color_mat[i,j] <- 15
      } else if(abs(theta - 1/5) <= 1e-3){
        color_mat[i,j] <- 16
      } else if(abs(theta - 2/5) <= 1e-3){
        color_mat[i,j] <- 17
      } else if(abs(theta - 3/5) <= 1e-3){
        color_mat[i,j] <- 18
      } else if(abs(theta - 4/5) <= 1e-3){
        color_mat[i,j] <- 19
      } else {
        color_mat[i,j] <- 20
      }
      
    } else {
      color_mat[i,j] <- 21
    }
  }
}

chromosome_value <- sapply(colnames(mat_theta), function(x){
  x <- strsplit(x = x, split = "_")[[1]][2]
  x <- strsplit(x = x, split = ":")[[1]][1]
  as.numeric(x)
})
chromosome_break_idx <- which(diff(chromosome_value) != 0)

clone_list <- lapply(1:6, function(x){
  which(SNU$clone == x)
})
clone_list[[7]] <- which(is.na(SNU$clone))
clone_len <- sapply(clone_list, length)
cell_idx <- unlist(clone_list)

color_mat <- color_mat[cell_idx,]

png("../../../../out/figures/Writeup14e/Writeup14e_SNU_atac_rawheatmap.png",
    height = 2000, width = 3000, units = "px", res = 300)
image(multiomicCCA:::.rotate(color_mat),
      breaks = seq(0.5, 21.5, by = 1), 
      col = col_vec)
for(i in 1:(length(clone_len)-1)){
  lines(x = c(-10, 10), y = rep(1-sum(clone_len[1:i])/sum(clone_len), 2),
        lwd = 2)
}
shift_val <- -1/(2*(ncol(color_mat)-1))
multiplier_val <- 2/(2*(ncol(color_mat)-1))
for(i in 1:(length(chromosome_break_idx))){
  lines(y = c(-10, 10), x = rep(shift_val+chromosome_break_idx[i]*multiplier_val, 2),
        lwd = 2)
}
graphics.off()

#########################


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

mat_theta <- pmin(mat_theta, quantile(mat_theta, probs = 0.95))
mat_rho <- pmin(mat_rho, quantile(mat_rho, probs = 0.95))


png("../../../../out/figures/Writeup14e/Writeup14e_SNU_atac_rawheatmap_theta.png",
    height = 2000, width = 3000, units = "px", res = 300)
image(multiomicCCA:::.rotate(mat_theta))
for(i in 1:(length(clone_len)-1)){
  lines(x = c(-10, 10), y = rep(1-sum(clone_len[1:i])/sum(clone_len), 2),
        lwd = 2)
}
shift_val <- -1/(2*(ncol(color_mat)-1))
multiplier_val <- 2/(2*(ncol(color_mat)-1))
for(i in 1:(length(chromosome_break_idx))){
  lines(y = c(-10, 10), x = rep(shift_val+chromosome_break_idx[i]*multiplier_val, 2),
        lwd = 2)
}
graphics.off()


png("../../../../out/figures/Writeup14e/Writeup14e_SNU_atac_rawheatmap_rho.png",
    height = 2000, width = 3000, units = "px", res = 300)
image(multiomicCCA:::.rotate(mat_rho))
for(i in 1:(length(clone_len)-1)){
  lines(x = c(-10, 10), y = rep(1-sum(clone_len[1:i])/sum(clone_len), 2),
        lwd = 2)
}
shift_val <- -1/(2*(ncol(color_mat)-1))
multiplier_val <- 2/(2*(ncol(color_mat)-1))
for(i in 1:(length(chromosome_break_idx))){
  lines(y = c(-10, 10), x = rep(shift_val+chromosome_break_idx[i]*multiplier_val, 2),
        lwd = 2)
}
graphics.off()

###################

set.seed(10)
svd_res <- irlba::irlba(mat_theta, nv = 50)
cna_k <- 7
dimred <-  multiomicCCA:::.mult_mat_vec(svd_res$u[,1:cna_k], svd_res$d[1:cna_k])
rownames(dimred) <- colnames(SNU)
colnames(dimred) <- paste0("pca_", 1:ncol(dimred))
set.seed(10)
umap_res <- Seurat::RunUMAP(dimred)
SNU[["theta"]] <- Seurat::CreateDimReducObject(umap_res@cell.embeddings, assay = "CNA")
plot1 <- Seurat::DimPlot(SNU, reduction = "theta",
                         group.by = "clone", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU: Theta"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_cna_theta.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")


set.seed(10)
svd_res <- irlba::irlba(mat_rho, nv = 50)
cna_k <- 7
dimred <-  multiomicCCA:::.mult_mat_vec(svd_res$u[,1:cna_k], svd_res$d[1:cna_k])
rownames(dimred) <- colnames(SNU)
colnames(dimred) <- paste0("pca_", 1:ncol(dimred))
set.seed(10)
umap_res <- Seurat::RunUMAP(dimred)
SNU[["rho"]] <- Seurat::CreateDimReducObject(umap_res@cell.embeddings, assay = "CNA")
plot1 <- Seurat::DimPlot(SNU, reduction = "rho",
                         group.by = "clone", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU: Rho"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_cna_rho.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

             
             