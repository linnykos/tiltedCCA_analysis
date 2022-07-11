rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_common <- multiSVD_obj$common_mat_1[cell_idx[order_vec],]

p <- ncol(rna_common)
predicted_list <- vector("list", length = p)
for(j in 1:p){
  if(j > 10 & j %% floor(p/10) == 0) {
    save(rna_common, predicted_list, 
         file = "../../../../out/Writeup14o/Writeup14o_10x_greenleaf_timeseries_tmp.RData")
  }
  print(j)
  
  genlasso_res <- genlasso::trendfilter(rna_common[,j], ord = 2)
  
  # https://stackoverflow.com/questions/34208564/how-to-hide-or-disable-in-function-printed-message
  invisible(capture.output(genlasso_cv <- genlasso::cv.trendfilter(genlasso_res)))
  res <- genlasso::predict.genlasso(genlasso_res, lambda = genlasso_cv$lambda.min)
  predicted_vec <- res$fit[,1]
  predicted_list[[j]] <- predicted_vec
}

save(greenleaf, multiSVD_obj, rna_common, predicted_list, 
     date_of_run, session_info,
     file = "../../../../out/Writeup14o/Writeup14o_10x_greenleaf_timeseries.RData")

# png(paste0("../../../../out/figures/Writeup14o/Writeup14o_10x_greenleaf_timeseries.png"),
#     height = 2000, width = 5000, units = "px", res = 300)
# par(mfrow = c(2,5), mar = c(4,4,4,0.5))
# for(j in 1:10){
#   plot(rna_common[,j], pch = 16, col = rgb(0.5, 0.5, 0.5, 0.2), main = j)
#   lines(predicted_mat[,j], col = "white", lwd = 4)
#   lines(predicted_mat[,j], col = 2, lwd = 2)
# }
# graphics.off()

