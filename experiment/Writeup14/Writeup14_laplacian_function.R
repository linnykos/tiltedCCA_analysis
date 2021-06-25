lap_eig <- function(mat, k_max = 200){
  res <- RSpectra::eigs(mat, k = k_max)
  res$values <- Re(res$values)
  res$vectors <- Re(res$vectors)
  inv_vectors <- MASS::ginv(res$vectors)
  
  approx_mat <- multiomicCCA:::.mult_mat_vec(res$vectors, res$values) %*% inv_vectors
  print(sqrt(sum((mat - approx_mat)^2))/sqrt(sum(mat^2)))
  
  tmp <- multiomicCCA:::.mult_mat_vec(res$vectors[,-1], res$values[-1])
  for(j in 1:ncol(tmp)){
    max_pos <- max(pmax(tmp[,j], 0))
    max_neg <- abs(min(pmin(tmp[,j],0)))
    if(max_neg > max_pos) tmp[,j] <- -tmp[,j]
  }
  
  tmp
}

compute_lap <- function(g, k_max = 100, normalize = T,
                        rowname_vec, colname_vec){
  n <- nrow(g)
  if(normalize){
    deg_vec <- sparseMatrixStats::rowSums2(g)
    invdeg_mat <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1/deg_vec)
    g <- invdeg_mat %*% g %*% invdeg_mat
  }
  
  deg_vec2 <- sparseMatrixStats::rowSums2(g)
  invdeg_mat2 <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1/deg_vec2)
  lap_mat <- invdeg_mat2 %*% g 
  
  basis <- lap_eig(lap_mat, k_max = k_max)
  rownames(basis) <- rowname_vec
  colnames(basis) <- colname_vec
  
  basis
}

compute_smooth_signal <- function(vec, basis){
  lm_res <- stats::lm(vec ~ basis) 
  list(pred_vec = lm_res$fitted.values, r_squared = summary(lm_res)$adj.r.squared)
}

create_plot <- function(seurat_obj, var_name, e_vec, c_vec, d_vec, 
                        e_res, c_res, d_res,
                        filename){
  gene_mat <- cbind(e_vec, c_vec, d_vec, e_res$pred_vec, c_res$pred_vec, d_res$pred_vec)
  rownames(gene_mat) <- colnames(seurat_obj)
  colnames(gene_mat) <- paste0("gene", 1:6)
  seurat_obj[["gene"]] <- Seurat::CreateDimReducObject(embedding = gene_mat, 
                                                       key = "gene", assay = "RNA")
  
  e_range <- range(c(e_vec, e_res$pred_vec))
  c_range <- range(c(c_vec, c_res$pred_vec))
  d_range <- range(c(d_vec, d_res$pred_vec))
  
  plot1 <- Seurat::FeaturePlot(seurat_obj, features = "gene_1", reduction = "everything")
  plot1 <- plot1 + ggplot2::labs(title = paste0(var_name, "\nRNA Denoised, Everything"), x = "Everything 1", y = "Everything 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = e_range)
  plot2 <- Seurat::FeaturePlot(seurat_obj, features = "gene_2", reduction = "common")
  plot2 <- plot2 + ggplot2::labs(title = paste0(var_name, "\nRNA Denoised, Common"), x = "Common 1", y = "Common 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = c_range)
  plot3 <- Seurat::FeaturePlot(seurat_obj, features = "gene_3", reduction = "distinct")
  plot3 <- plot3 + ggplot2::labs(title = paste0(var_name, "\nRNA Denoised, Distinct"), x = "Distinct 1", y = "Distinct 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = d_range)

  plot4 <- Seurat::FeaturePlot(seurat_obj, features = "gene_4", reduction = "everything")
  plot4 <- plot4 + ggplot2::labs(title = paste0("Smoothed Everything\nVar: ", round(stats::var(e_res$pred_vec), 2), ", R2: ", round(e_res$r_squared, 2)), 
                                 x = "Everything 1", y = "Everything 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = e_range)
  plot5 <- Seurat::FeaturePlot(seurat_obj, features = "gene_5", reduction = "common")
  plot5 <- plot5 + ggplot2::labs(title = paste0("Smoothed Common\nVar: ", round(stats::var(c_res$pred_vec), 2), ", R2: ", round(c_res$r_squared, 2)), 
                                 x = "Common 1", y = "Common 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = c_range)
  plot6 <- Seurat::FeaturePlot(seurat_obj, features = "gene_6", reduction = "distinct")
  plot6 <- plot6 + ggplot2::labs(title = paste0("Smoothed Distinct\nVar: ", round(stats::var(d_res$pred_vec), 2), ", R2: ", round(d_res$r_squared, 2)), 
                                 x = "Distinct 1", y = "Distinct 2") + 
    ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = d_range)
  
  p <- cowplot::plot_grid(plot1, plot2, plot3, plot4, plot5, plot6)
  cowplot::save_plot(filename = filename, p, ncol = 3, nrow = 2, base_asp = 1.2, device = "png")
  
  invisible()
}
