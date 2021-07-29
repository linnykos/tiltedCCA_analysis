rm(list=ls())
load("../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")

nn <- 30
set.seed(10)
combined_g <- multiomicCCA::combine_frnn(dcca_res, g_1 = rna_frnn$c_g,
                                         g_2 = atac_frnn$c_g, nn = nn, 
                                         verbose = 2)

tmp <- myeloid2[["distinct"]]@cell.embeddings
# idx1 <- which(myeloid2@meta.data$celltype == "B16_dICB")
idx1 <- which(tmp[,1] <= -3.5)
idx2 <- which(tmp[,2] >= -2.5)
idx3 <- which(tmp[,2] <= -1)
start_idx <- intersect(idx1, intersect(idx2, idx3))
table(myeloid2@meta.data$celltype[start_idx])

.clisi_cell <- function(g, idx, tol = 1e-3){
  target_bg <- length(idx)/nrow(g)
  
  vec <- sapply(idx, function(i){
    neigh <- multiomicCCA:::.nonzero_col(g, i, bool_value = F)
    len <- length(neigh)
    if(len == 0){
      return(0)
    }
    in_len <- length(which(neigh %in% idx))
    
    max((in_len/len - target_bg + tol)/(1-target_bg+tol), 0)
  })
  
  median(vec)
}

common_score <- .clisi_cell(combined_g, start_idx)
distinct_score1 <- .clisi_cell(rna_frnn$d_g, start_idx)
distinct_score2 <- .clisi_cell(atac_frnn$d_g, start_idx)
base_vec <- c(common_score, distinct_score1, distinct_score2)
tol <- 0.02
base_vec

iter <- 1
tries <- 10
idx <- start_idx
while(TRUE){
  print(paste0("Iteration ", iter, ": Length of ", length(idx)))
  neigh_list <- unlist(lapply(idx, function(i){
    multiomicCCA:::.nonzero_col(combined_g, i, bool_value = F)
  }))
  neigh_list <- setdiff(neigh_list, idx)
  if(length(neigh_list) == 0) break()
  candidates <- sort(table(neigh_list), decreasing = T)
  try_idx <- 1
  
  while(try_idx <= tries){
    bool_continue <- FALSE
    tmp_idx <- c(idx, as.numeric(names(candidates[try_idx])))
    common_score <- .clisi_cell(combined_g, tmp_idx)
    distinct_score1 <- .clisi_cell(rna_frnn$d_g, tmp_idx)
    distinct_score2 <- .clisi_cell(atac_frnn$d_g, tmp_idx)
    
    bool1 <- common_score <= base_vec[1]+tol
    bool2 <- distinct_score1 >= 2/3 #base_vec[2]-tol
    bool3 <- distinct_score2 <= base_vec[3]+tol
    
    if(bool1 & bool2 & bool3){
      idx <- unique(tmp_idx); bool_continue <- TRUE; break()
    } else {
      try_idx <- try_idx + 1
    }
  }
  
  if(!bool_continue) break()
  
  print(round(c(common_score, distinct_score1, distinct_score2), 2))
  iter <- iter + 1
}
