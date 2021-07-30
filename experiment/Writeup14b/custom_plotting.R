.plot_umap <- function(layout, full_col_vec, idx, main){
  tmp_col_vec <- full_col_vec; tmp_col_vec[-idx] <- "gray"
  size_vec <- rep(0.5, n); size_vec[idx] <- 2
  plot(layout[-idx,1], layout[-idx,2], 
       xlim = range(layout[,1]), ylim = range(layout[,2]),
       asp = T, col = tmp_col_vec[-idx], 
       cex = size_vec[-idx], pch = 16, xlab = "UMAP 1", ylab = "UMAP 2",
       main = main)
  points(layout[idx,1], layout[idx,2],
         col = tmp_col_vec[idx], cex = size_vec[idx], pch = 16)
  invisible()
}

.custom_plot1 <- function(myeloid2, full_col_vec, res){
  par(mfcol = c(2,3), mar = c(4,4,4,0.1))
  layout <- myeloid2[["combined"]]@cell.embeddings
  .plot_umap(layout, full_col_vec, res$baseline_idx, 
             main = paste0("Common combined,\nInitial (", round(res$baseline_scores[1], 2), ")"))
  .plot_umap(layout, full_col_vec, res$idx, 
             main = paste0("Common combined,\nFinal (", round(res$scores[1], 2), ")"))
  
  layout <- myeloid2[["distinct"]]@cell.embeddings
  .plot_umap(layout, full_col_vec, res$baseline_idx, 
             main = paste0("RNA distinct,\nInitial (", round(res$baseline_scores[2], 2), ")"))
  .plot_umap(layout, full_col_vec, res$idx, 
             main = paste0("RNA distinct,\nFinal (", round(res$scores[2], 2), ")"))
  
  layout <- myeloid2[["distinct2"]]@cell.embeddings
  .plot_umap(layout, full_col_vec, res$baseline_idx,
             main = paste0("ATAC distinct,\nInitial (", round(res$baseline_scores[3], 2), ")"))
  .plot_umap(layout, full_col_vec, res$idx,
             main = paste0("ATAC distinct,\nFinal (", round(res$scores[3], 2), ")"))
  
  invisible()
}

.custom_plot2 <- function(myeloid2, full_col_vec, res){
  par(mfcol = c(2,3), mar = c(4,4,4,0.1))
  layout <- myeloid2[["everything"]]@cell.embeddings
  .plot_umap(layout, full_col_vec, res$baseline_idx, "RNA everything,\nInitial")
  .plot_umap(layout, full_col_vec, res$idx, "RNA everything,\nFinal")
  
  layout <- myeloid2[["everything2"]]@cell.embeddings
  .plot_umap(layout, full_col_vec, res$baseline_idx, "ATAC everything,\nInitial")
  .plot_umap(layout, full_col_vec, res$idx, "ATAC everything,\nFinal")
  
  layout <- myeloid2[["both"]]@cell.embeddings
  .plot_umap(layout, full_col_vec, res$baseline_idx, "Both everything,\nInitial")
  .plot_umap(layout, full_col_vec, res$idx, "Both everything,\nFinal")
  
  invisible()
}
