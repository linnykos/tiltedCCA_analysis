rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_greenleaf_timeseries.RData")

p <- ncol(rna_common)
sd_vec <- sapply(1:p, function(j){
  stats::sd(rna_predicted_mat[,j] - rna_common[,j])
})
sd_mean_vec <- sapply(1:p, function(j){
  stats::sd(rna_common[,j])
})

n <- nrow(rna_predicted_mat)
rna_category_vec <- rep(NA, p)
rna_category_vec[which(sd_vec/sd_mean_vec >= 0.9)] <- "flat"
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

atac_category_vec <- rep(NA, p)
atac_category_vec[which(sd_vec/sd_mean_vec >= 0.9)] <- "flat"
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


