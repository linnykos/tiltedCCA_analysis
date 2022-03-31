rm(list=ls())
load("../../../../out/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_tcca.RData")
library(Seurat)
source("../Writeup14l/slingshot.R")
source("nonparametric_curve.R")

set.seed(10)
greenleaf <- Seurat::FindMultiModalNeighbors(greenleaf, reduction.list = list("pca", "lsi"), 
                                             dims.list = list(1:50, 2:50))
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf,
                                  graph.name = "wsnn", algorithm = 3, 
                                  resolution = 2)
lineage_order <- as.character(c(11,15,8,3,14,0,1,2,6,12))

####

multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
dimred_1 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, 
                                     apply_postDimred = T, what = "common_mat")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
dimred_2 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, what = "common_mat")
dimred <- cbind(dimred_1, dimred_2)
svd_tmp <- irlba::irlba(dimred, nv = 25)
dimred <- tiltedCCA:::.mult_mat_vec(svd_tmp$u[,1:10], svd_tmp$d[1:10])

cluster_vec <- greenleaf$seurat_clusters[which(greenleaf$seurat_clusters %in% lineage_order)]
dimred <- dimred[which(greenleaf$seurat_clusters %in% lineage_order),]

initial_fit <- .initial_curve_fit(cluster_vec = cluster_vec,
                                  dimred = dimred,
                                  lineage_order = lineage_order)
pseudotime_vec <- .extract_pseudotime(dimred = dimred,
                                      initial_fit = initial_fit,
                                      stretch = 2)
n <- ncol(greenleaf)
pseudotime_full <- rep(NA, n)
pseudotime_full[which(greenleaf$seurat_clusters %in% lineage_order)] <- pseudotime_vec
pseudotime_full2 <- rep(NA, n)
pseudotime_full2[which(greenleaf$seurat_clusters %in% lineage_order)] <- rank(pseudotime_vec)

greenleaf$pseudotime <- pseudotime_full
greenleaf$pseudotime_rank <- pseudotime_full2

####

rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(rna = rna_common[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})

greenleaf$alignment <- alignment_vec

####

rna_common <- multiSVD_obj$common_mat_1   
atac_everything <- multiSVD_obj$common_mat_2 + multiSVD_obj$distinct_mat_2

gene_names <- sapply(colnames(atac_everything), function(x){
  tmp <- strsplit(x, split="-")[[1]]
  paste0(tmp[2:length(tmp)], collapse="-")
})
names(gene_names) <- NULL
all(gene_names %in% colnames(rna_common))
rna_common <- rna_common[,gene_names]
all(colnames(rna_common) == colnames(gene_names))

cell_idx <- which(greenleaf$seurat_clusters %in% lineage_order)
rna_common <- rna_common[cell_idx,]
atac_everything <- atac_everything[cell_idx,]
rna_common <- scale(rna_common)
atac_everything <- scale(atac_everything)

####

set.seed(10)
celltype_vec <- greenleaf$celltype[cell_idx]
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
pseudotime_vec <- pseudotime_vec + runif(length(pseudotime_vec), min = 0, max = 0.001)
order_vec <- order(pseudotime_vec, decreasing = F)
alignment_vec <- greenleaf$alignment[cell_idx]

rna_common2 <- rna_common[order_vec,]
atac_everything2 <- atac_everything[order_vec,]
celltype_vec <- celltype_vec[order_vec]
alignment_vec <- alignment_vec[order_vec]

####

p <- ncol(rna_common2)
n <- nrow(rna_common2)
rna_fit_list <- lapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  .shape_constrained_fit_search(rna_common2[,j],
                                grid_vec = round(seq(0,n,length.out=7)[-c(1,7)]),
                                metric = 2)
})
table(sapply(rna_fit_list, function(x){x$type}))
quantile(sapply(rna_fit_list, function(x){x$quality}))

atac_fit_list <- lapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  .shape_constrained_fit_search(atac_everything2[,j],
                                grid_vec = round(seq(0,n,length.out=7)[-c(1,7)]),
                                metric = 2)
})
quantile(sapply(atac_fit_list, function(x){x$quality}))

table(sapply(atac_fit_list, function(x){x$type}),
      sapply(rna_fit_list, function(x){x$type}))

###########################

tab_mat <- table(sapply(atac_fit_list, function(x){x$type}),
                 sapply(rna_fit_list, function(x){x$type}))

n <- nrow(rna_common2)
num_considered <- 100
for(atac_setting in rownames(tab_mat)){
  print(paste0("ATAC: ", atac_setting))
  atac_idx <- which(sapply(atac_fit_list, function(x){x$type}) == atac_setting)
  atac_subidx <- atac_idx[order(sapply(atac_fit_list[atac_idx], function(x){x$quality}), decreasing = T)]
  
  for(rna_setting in colnames(tab_mat)){
    print(paste0("RNA: ", rna_setting))
    rna_idx <- which(sapply(rna_fit_list, function(x){x$type}) == rna_setting)
    rna_subidx <- rna_idx[order(sapply(rna_fit_list[rna_idx], function(x){x$quality}), decreasing = T)]
    
    set.seed(10)
    comb_idx <- intersect(atac_subidx[1:num_considered], rna_subidx[1:num_considered])
    if(length(comb_idx) < 15){
      comb_idx <- intersect(atac_subidx, rna_subidx)
    }
    
    # first plot
    col_vec <- grDevices::colorRampPalette(c(3, 'lightgray', 2))(n)
    png(paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivityRelation_atac-",
               atac_setting, "_rna-", rna_setting, ".png"),
        height = 2000, width = 3000, units = "px", res = 300)
    par(mfrow = c(3,5), mar = c(4,4,4,0.5))
    for(j in comb_idx[1:min(15,length(comb_idx))]){
      shuff_idx <- sample(1:n)
      plot(atac_everything2[shuff_idx,j], rna_common2[shuff_idx,j],
           main = paste0("Correlation: ", round(stats::cor(atac_everything2[,j], rna_common2[,j]), 2)),
           xlab = paste0("Gene activity for ", colnames(rna_common2)[j]),
           ylab = paste0("Common RNA for ", colnames(rna_common2)[j]),
           col = col_vec[shuff_idx], pch = 16, cex = 0.5
      )
    }
    graphics.off()
    
    # second plot
    col_vec <- grDevices::colorRampPalette(c(3, 'lightgray', 2))(n)
    png(paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivityRelation_atac-",
               atac_setting, "_rna-", rna_setting, "_fitted.png"),
        height = 2000, width = 3000, units = "px", res = 300)
    par(mfrow = c(3,5), mar = c(4,4,4,0.5))
    for(j in comb_idx[1:min(15,length(comb_idx))]){
      shuff_idx <- sample(1:n)
      plot(atac_fit_list[[j]]$fit[shuff_idx], rna_fit_list[[j]]$fit[shuff_idx],
           main = paste0("ATAC: ", atac_setting, ",\nRNA: ", rna_setting),
           xlab = paste0("Gene activity for ", colnames(rna_common2)[j]),
           ylab = paste0("Common RNA for ", colnames(rna_common2)[j]),
           col = col_vec[shuff_idx], pch = 16, cex = 2
      )
    }
    graphics.off()
    
    # third plot
    col_pal <- grDevices::colorRampPalette(c('lightgrey', 'blue'))(100)
    seq_vec <- seq(min(alignment_vec), max(alignment_vec), length.out = 100)
    col_vec <- sapply(alignment_vec, function(val){
      col_pal[which.min(abs(val - seq_vec))]
    })
    png(paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivityRelation_atac-",
               atac_setting, "_rna-", rna_setting, "_alignment.png"),
        height = 2000, width = 3000, units = "px", res = 300)
    par(mfrow = c(3,5), mar = c(4,4,4,0.5))
    for(j in comb_idx[1:min(15,length(comb_idx))]){
      plot(atac_everything2[,j], rna_common2[,j],
           xlab = paste0("Gene activity for ", colnames(rna_common2)[j]),
           ylab = paste0("Common RNA for ", colnames(rna_common2)[j]),
           col = col_vec, pch = 16, cex = 0.5
      )
    }
    graphics.off()
    
    # fourth plot
    png(paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivityRelation_atac-",
               atac_setting, "_rna-", rna_setting, "_vsTime.png"),
        height = 2000, width = 3000, units = "px", res = 300)
    par(mfrow = c(3,5), mar = c(4,4,4,0.5))
    for(j in comb_idx[1:min(15,length(comb_idx))]){
      plot(NA, 
           xlim = c(1,n),
           ylim = range(c(atac_everything2[,j], rna_common2[,j])),
           xlab = "Rank (via Slingshot)",
           ylab = paste0("Expression for ", colnames(rna_common2)[j])
      )
      
      x_vec <- c(1:n, 1:n)
      y_vec <- c(rna_common2[,j], atac_everything2[,j])
      col_vec <- c(rep(rgb(0.5, 0.5, 0.5, 0.5), n),
                   rep(rgb(0.75, 0, 0, 0.5), n))
      shuff_idx <- sample(1:length(x_vec))
      points(x = x_vec[shuff_idx], y = y_vec[shuff_idx], 
             col = col_vec[shuff_idx], pch = 16, cex = 0.5)
      lines(x = 1:n, y = rna_fit_list[[j]]$fit,
            col = "white", lwd = 5)
      lines(x = 1:n, y = rna_fit_list[[j]]$fit,
            col = "black", lwd = 3)
      lines(x = 1:n, y = atac_fit_list[[j]]$fit,
            col = "white", lwd = 5)
      lines(x = 1:n, y = atac_fit_list[[j]]$fit,
            col = 2, lwd = 3)
    }
    graphics.off()
  }
}

#####################################################

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
idx <- which(colnames(rna_common2) %in% c(s.genes, g2m.genes))

table(sapply(atac_fit_list[idx], function(x){x$type}),
      sapply(rna_fit_list[idx], function(x){x$type}))

col_vec <- grDevices::colorRampPalette(c(3, 'lightgray', 2))(n)
png(paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_geneactivityRelation_cellcyle_fitted.png"),
    height = 2500, width = 3000, units = "px", res = 300)
par(mfrow = c(5,6), mar = c(4,4,4,0.5))
for(j in idx[1:30]){
  shuff_idx <- sample(1:n)
  plot(atac_fit_list[[j]]$fit[shuff_idx], rna_fit_list[[j]]$fit[shuff_idx],
       main = paste0("ATAC: ", atac_fit_list[[j]]$type, ",\nRNA: ", rna_fit_list[[j]]$type),
       xlab = paste0("Gene activity for ", colnames(rna_common2)[j]),
       ylab = paste0("Common RNA for ", colnames(rna_common2)[j]),
       col = col_vec[shuff_idx], pch = 16, cex = 2
  )
}
graphics.off()
