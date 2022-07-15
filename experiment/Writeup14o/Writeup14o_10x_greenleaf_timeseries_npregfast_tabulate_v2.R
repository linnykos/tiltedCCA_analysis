rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

# load("../../../../out/main/10x_greenleaf_timeseries.RData")

# save(atac_pred, atac_predicted_mat, 
#      rna_common, rna_predicted_mat, 
#      file = "../../../../out/tmp.RData")
load("../../out/tmp.RData")

colnames(atac_predicted_mat) <- colnames(atac_pred)
colnames(rna_predicted_mat) <- colnames(rna_common)

# j <- which(colnames(rna_common) == "IL1A")
# j <- which(colnames(rna_common) == "CCL5")
# j <- which(colnames(rna_common) == "PLA2G7")
j <- which(colnames(rna_common) == "S100B")

n <- nrow(rna_common)
x_vec <- c(1:n)/n
tmp_df <- data.frame(y = atac_predicted_mat[,j], x = x_vec)
lm_res <- stats::lm(y ~ poly(x,5,raw = T), data = tmp_df)
coef_vec <- stats::coef(lm_res)

plot(x_vec, atac_predicted_mat[,j])
points(x_vec, lm_res$fitted.values, col = 2, pch = 16)
# points(x_vec, coef_vec[1] + coef_vec[2]*x_vec + coef_vec[3]*x_vec^2 + coef_vec[4]*x_vec^3, col = 3, pch = 16)

####

assess_curve <- function(vec, 
                         range_vec = c(0,1),
                         tol = 1e-4){
  n <- length(vec)
  x_vec <- c(1:n)/n
  tmp_df <- data.frame(y = vec, x = x_vec)
  lm_res <- stats::lm(y ~ poly(x,5,raw = T), data = tmp_df)
  coef_vec <- stats::coef(lm_res)
  
  coef_vec_gradient <- coef_vec[-1]*c(1,2,3,4,5)
  critical <- polyroot(coef_vec_gradient)
  critical <- critical[any(abs(Im(critical)) <= tol)]
  critical <- Re(critical)
  
  labeling <- sapply(critical, function(val){
    res <- coef_vec_gradient %*% c(1,val,val^2,val^3,val^4)
    if(abs(res) <= tol){
      res2 <- (coef_vec_gradient[-1]*c(1,2,3,4)) %*% c(1,val,val^2,val^3)
      if(res2 > 0) return("loc_min") else return("loc_max")
    } else {
      return("saddle")
    }
  })
  
  x_vec <- c(range_vec, critical)
  type_vec <- c("boundary", "boundary", unlist(labeling))
  df <- data.frame(x = x_vec, type = type_vec)
  df$y <- sapply(df$x, function(x){
    coef_vec %*% c(1,x,x^2,x^3,x^4,x^5)
  })
  df <- df[,c("x", "y", "type")]
  df <- df[order(df$x),]
  
  df
}

assess_curve(atac_predicted_mat[,"MT-CO2"])
assess_curve(atac_predicted_mat[,"AC011029.1"])
assess_curve(atac_predicted_mat[,"S100B"])
assess_curve(atac_predicted_mat[,"RHOJ"])
assess_curve(atac_predicted_mat[,"IL1A"])
assess_curve(atac_predicted_mat[,"PLA2G7"])
assess_curve(atac_predicted_mat[,"DLX5"])
assess_curve(atac_predicted_mat[,"CCL4L2"])
assess_curve(atac_predicted_mat[,"CADM1"])
assess_curve(atac_predicted_mat[,"AC005498.1"])
assess_curve(atac_predicted_mat[,"CNTN4"])
assess_curve(atac_predicted_mat[,"PNOC"])
assess_curve(rna_predicted_mat[,"CSF2RA"])
assess_curve(atac_predicted_mat[,"SCGN"])
assess_curve(rna_predicted_mat[,"KIF18B"])
assess_curve(atac_predicted_mat[,"PTPRG-AS1"])
assess_curve(atac_predicted_mat[,"ARHGAP24"])
assess_curve(atac_predicted_mat[,"CKS2"])


