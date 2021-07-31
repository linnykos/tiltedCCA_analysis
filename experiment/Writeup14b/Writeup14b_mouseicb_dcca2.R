rm(list = ls())

library(Seurat)
library(multiomicCCA)

load("../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)
mat_1_denoised <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
mat_2_denoised <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2

rm(list = "dcca_decomp"); gc(T)
c_eig <- myeloid2[["clap"]]@cell.embeddings
d_eig <- myeloid2[["dlap"]]@cell.embeddings
e_eig <- myeloid2[["elap"]]@cell.embeddings

print("Starting RNA smooth")
p1 <- ncol(mat_1_denoised)
gene_smoothed <- lapply(1:p1, function(j){
  if(j %% floor(p1/10) == 0) cat('*')
  
  c_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,j], c_eig)
  d_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,j], d_eig)
  e_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,j], e_eig)
  
  list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
       d_variance = d_res$variance, d_r2 = d_res$r_squared,
       e_variance = e_res$variance, e_r2 = e_res$r_squared)
})

save(date_of_run, session_info, dcca_res, myeloid2, rna_frnn, atac_frnn, combined_g,
     gene_smoothed, 
     file = "../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")

#########

c_eig2 <- myeloid2[["clap2"]]@cell.embeddings
d_eig2 <- myeloid2[["dlap2"]]@cell.embeddings
e_eig2 <- myeloid2[["elap2"]]@cell.embeddings

print("Starting ATAC smooth")
p2 <- ncol(mat_2_denoised)
sd_vec <- matrixStats::colSds(mat_2_denoised)
peak_order <- order(sd_vec, decreasing = T)
num_breaks <- 100
atac_smoothed <- vector("list", p2)
break_vec <- round(seq(1, p2+1, length.out = num_breaks+1))
for(x in 1:(length(break_vec)-1)){
  print(paste0("On subset ", x))
  tmp_idx <- peak_order[break_vec[x]:(break_vec[x+1]-1)]
  
  tmp_list <- lapply(tmp_idx, function(j){
    c_res <- multiomicCCA::compute_smooth_signal(mat_2_denoised[,j], c_eig2)
    d_res <- multiomicCCA::compute_smooth_signal(mat_2_denoised[,j], d_eig2)
    e_res <- multiomicCCA::compute_smooth_signal(mat_2_denoised[,j], e_eig2)
    
    list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
         d_variance = d_res$variance, d_r2 = d_res$r_squared,
         e_variance = e_res$variance, e_r2 = e_res$r_squared)
  })
  
  atac_smoothed[tmp_idx] <- tmp_list
  
  save(date_of_run, session_info, dcca_res, myeloid2, rna_frnn, atac_frnn, combined_g,
       gene_smoothed, atac_smoothed,
       file = "../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")
}


save(date_of_run, session_info, dcca_res, myeloid2, rna_frnn, atac_frnn, combined_g,
     gene_smoothed, atac_smoothed,
     file = "../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")
