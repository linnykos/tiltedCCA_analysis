rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/main/10x_greenleaf_timeseries.RData")

p <- ncol(rna_common)
sd_vec <- sapply(1:p, function(j){
  stats::sd(rna_predicted_mat[,j] - rna_common[,j])
})
sd_mean_vec <- sapply(1:p, function(j){
  stats::sd(rna_common[,j])
})
flat_idx <- intersect(which(sd_vec/sd_mean_vec >= 0.9),
                      which(sapply(1:p, function(j){diff(range(rna_predicted_mat[,j])) <= sd_vec[j]})))

n <- nrow(rna_predicted_mat)
rna_category_vec <- rep(NA, p)
rna_category_vec[flat_idx] <- "flat"
rna_category_vec2 <- sapply(1:p, function(j){
  if(j > 10 && j %% floor(p/10) == 0) cat('*')
  if(!is.na(rna_category_vec[j]) && rna_category_vec[j] == "flat") return("flat")
  
  vec <- rna_predicted_mat[,j]
  res_list <- vector("list", 4)
  names(res_list) <- c("increasing", "decreasing", "upward_bump", "downward_bump")
  
  res_list[["increasing"]] <- UniIsoRegression::reg_1d(
    vec,
    w_vec = rep(1, n),
    metric = 2,
    unimodal = F,
    decreasing = F)
  
  res_list[["decreasing"]] <- UniIsoRegression::reg_1d(
    vec,
    w_vec = rep(1, n),
    metric = 2,
    unimodal = F,
    decreasing = T)
  
  res_list[["upward_bump"]] <- UniIsoRegression::reg_1d(
    vec,
    w_vec = rep(1, n),
    metric = 2,
    unimodal = T,
    decreasing = F)
  res_list[["downward_bump"]] <- -1*UniIsoRegression::reg_1d(
    -vec,
    w_vec = rep(1, n),
    metric = 2,
    unimodal = T,
    decreasing = F)
  sd_curves <- sapply(res_list, function(vec2){stats::sd(vec - vec2)})
  if(which.min(sd_curves) == 1) return("increasing")
  if(which.min(sd_curves) == 2) return("decreasing")
  
  threshold <- 1
  if(which.min(sd_curves) == 3){
    max_val <- max(res_list[["upward_bump"]])
    max_idx <- which.max(res_list[["upward_bump"]])
    min_val_1 <- min(res_list[["upward_bump"]][1:max_idx])
    min_val_2 <- min(res_list[["upward_bump"]][max_idx:n])
    
    if(max_val - max(c(min_val_1, min_val_2)) >= threshold*sd_vec[j]){
      return("upward_bump")
    }
  } else if(which.min(sd_curves) == 4){
    min_val <- min(res_list[["downward_bump"]])
    min_idx <- which.min(res_list[["downward_bump"]])
    max_val_1 <- max(res_list[["downward_bump"]][1:min_idx])
    max_val_2 <- max(res_list[["downward_bump"]][min_idx:n])
    
    if(min(c(max_val_1, max_val_2)) - min_val >= threshold*sd_vec[j]){
      return("downward_bump")
    }
  }
  
  if(sd_curves[1] <= sd_curves[2]) return("increasing")
  if(sd_curves[1] >= sd_curves[2]) return("decreasing")
})

######################

sd_vec <- sapply(1:p, function(j){
  stats::sd(atac_predicted_mat[,j] - atac_pred[,j])
})
sd_mean_vec <- sapply(1:p, function(j){
  stats::sd(atac_pred[,j])
})
flat_idx <- intersect(which(sd_vec/sd_mean_vec >= 0.9),
                      which(sapply(1:p, function(j){diff(range(atac_predicted_mat[,j])) <= sd_vec[j]})))

atac_category_vec <- rep(NA, p)
atac_category_vec[flat_idx] <- "flat"
atac_category_vec2 <- sapply(1:p, function(j){
  if(j > 10 && j %% floor(p/10) == 0) cat('*')
  if(!is.na(atac_category_vec[j]) && atac_category_vec[j] == "flat") return("flat")
  
  vec <- atac_predicted_mat[,j]
  res_list <- vector("list", 4)
  names(res_list) <- c("increasing", "decreasing", "upward_bump", "downward_bump")
  
  res_list[["increasing"]] <- UniIsoRegression::reg_1d(
    vec,
    w_vec = rep(1, n),
    metric = 2,
    unimodal = F,
    decreasing = F)
  
  res_list[["decreasing"]] <- UniIsoRegression::reg_1d(
    vec,
    w_vec = rep(1, n),
    metric = 2,
    unimodal = F,
    decreasing = T)
  
  res_list[["upward_bump"]] <- UniIsoRegression::reg_1d(
    vec,
    w_vec = rep(1, n),
    metric = 2,
    unimodal = T,
    decreasing = F)
  res_list[["downward_bump"]] <- -1*UniIsoRegression::reg_1d(
    -vec,
    w_vec = rep(1, n),
    metric = 2,
    unimodal = T,
    decreasing = F)
  sd_curves <- sapply(res_list, function(vec2){stats::sd(vec - vec2)})
  if(which.min(sd_curves) == 1) return("increasing")
  if(which.min(sd_curves) == 2) return("decreasing")
  
  threshold <- 1
  if(which.min(sd_curves) == 3){
    max_val <- max(res_list[["upward_bump"]])
    max_idx <- which.max(res_list[["upward_bump"]])
    min_val_1 <- min(res_list[["upward_bump"]][1:max_idx])
    min_val_2 <- min(res_list[["upward_bump"]][max_idx:n])
    
    if(max_val - max(c(min_val_1, min_val_2)) >= threshold*sd_vec[j]){
      return("upward_bump")
    }
  } else if(which.min(sd_curves) == 4){
    min_val <- min(res_list[["downward_bump"]])
    min_idx <- which.min(res_list[["downward_bump"]])
    max_val_1 <- max(res_list[["downward_bump"]][1:min_idx])
    max_val_2 <- max(res_list[["downward_bump"]][min_idx:n])
    
    if(min(c(max_val_1, max_val_2)) - min_val >= threshold*sd_vec[j]){
      return("downward_bump")
    }
  }
  
  if(sd_curves[1] <= sd_curves[2]) return("increasing")
  if(sd_curves[1] >= sd_curves[2]) return("decreasing")
})

table(atac_category_vec2, rna_category_vec2)

#####################

uniq_type <- sort(unique(rna_category_vec2))
n <- nrow(rna_common)

for(i in uniq_type){
  print(paste0("RNA: ", i))
  for(j in uniq_type){
    print(paste0("ATAC: ", j))
    
    idx <- intersect(which(rna_category_vec2 == i), which(atac_category_vec2 == j))
    set.seed(10)
    if(length(idx) > 16) idx <- sample(idx, 16)
    
    png(paste0("../../../../out/figures/Writeup14o/Writeup14o_10x_greenleaf_timeseries_atac-",
               j, "_rna-", i, ".png"),
        height = 2000, width = 2000, units = "px", res = 300)
    par(mfrow = c(4,4), mar = c(4,4,4,0.5))
    for(k in idx){
      ylim <- quantile(c(rna_common[,k], atac_pred[,k]), probs = c(0.05, 0.95))
      plot(x = 1:n, y = rna_common[,k], col = rgb(0.5, 0.5, 0.5, 0.2), ylim = ylim,
           main = colnames(rna_common)[k], xlab = "Slingshot rank", ylab = "Value")
      points(x = 1:n, y = atac_pred[,k], col = rgb(0.8, 0, 0, 0.2))
      
      lines(x = 1:n, y = rna_predicted_mat[,k], col = "white", lwd = 4)
      lines(x = 1:n, y = atac_predicted_mat[,k], col = "white", lwd = 4)
      lines(x = 1:n, y = rna_predicted_mat[,k], col = 1, lwd = 2.5)
      lines(x = 1:n, y = atac_predicted_mat[,k], col = 2, lwd = 3)
    }
    graphics.off()
    
  }
}




