rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/Writeup14d/Writeup14d_10x_mouseembryo_dcca.RData")

source("diffusion_distance.R")
source("frnn_alt.R")
source("changepoint.R")
source("plotting.R")
source("trajectory.R")

clustering <- as.character(mbrain2@meta.data$wsnn_res.2)
selected_clusters <- c("17", "16", "0", "4", "2", "1", "5", "7", "12", "20")

set.seed(10)
res <- form_snn_graph(dcca_res, k = 5, 
                      clustering = clustering,
                      selected_clusters = selected_clusters)
P <- form_transition(res$snn,
                     lazy_param = 0.85,
                     teleport_param = 0.99)
diffusion_res <- extract_eigen(P, check = T)
n <- nrow(P)
dist_mat <- diffusion_distance(diffusion_res$eigenvalues, 
                               diffusion_res$right_vector,
                               time_vec = seq(1,5,by=1))
rownames(dist_mat) <- rownames(P)
colnames(dist_mat) <- colnames(P)

radial_start <- which(rownames(res$adj_mat) == "17")
cortical_end <- which.max(dist_mat[radial_start,])

# compute janky pseudotime
forward_rank <- rank(dist_mat[radial_start,])
backward_rank <- nrow(res$adj_mat) - rank(dist_mat[cortical_end,])
avg_rank <- apply(cbind(forward_rank, backward_rank), 1, mean)
avg_rank <- rank(avg_rank, ties.method = "first")

selected_clusters <- colnames(dist_mat)
selected_clusters_ordering <- avg_rank

prin_df <- compute_principal_ordering(dcca_res,
                                      clustering = clustering, 
                                      selected_clusters = selected_clusters,
                                      selected_clusters_ordering = selected_clusters_ordering)

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)
common_changepoint <- detect_changepoints(obj = dcca_decomp$common_mat_1,
                                          prin_df = prin_df,
                                          threshold = 0.01/ncol(dcca_decomp$common_mat_1))
distinct_changepoint <- detect_changepoints(obj = dcca_decomp$distinct_mat_1,
                                          prin_df = prin_df,
                                          threshold = 0.01/ncol(dcca_decomp$distinct_mat_1))

common_pval <- sapply(common_changepoint, function(x){x$pval})
quantile(common_pval)
distinct_pval <- sapply(distinct_changepoint, function(x){x$pval})
quantile(distinct_pval)

intersection_genes <- intersect(names(common_changepoint), names(distinct_changepoint))

# try removing all intersection based on height
common_changepoint2 <- common_changepoint
distinct_changepoint2 <- distinct_changepoint
for(gene in intersection_genes){
  idx1 <- which(names(common_changepoint2) == gene)
  idx2 <- which(names(distinct_changepoint2) == gene)
  if(common_changepoint2[[idx1]]$diff > distinct_changepoint2[[idx2]]$diff){
    distinct_changepoint2[[idx2]] <- NA
  } else {
    common_changepoint2[[idx1]] <- NA
  }
}
common_changepoint2 <- common_changepoint2[sapply(common_changepoint2, function(x){!all(is.na(x))})]
distinct_changepoint2 <- distinct_changepoint2[sapply(distinct_changepoint2, function(x){!all(is.na(x))})]

# try removing all intersection based on length
common_changepoint3 <- common_changepoint
distinct_changepoint3 <- distinct_changepoint
for(gene in intersection_genes){
  idx1 <- which(names(common_changepoint3) == gene)
  idx2 <- which(names(distinct_changepoint3) == gene)
  len1 <- common_changepoint3[[idx1]]$end - common_changepoint3[[idx1]]$start
  len2 <- distinct_changepoint3[[idx2]]$end - distinct_changepoint3[[idx2]]$start
  
  if(len2 > len1){
    distinct_changepoint3[[idx2]] <- NA
  } else {
    common_changepoint3[[idx1]] <- NA
  }
}
common_changepoint3 <- common_changepoint3[sapply(common_changepoint3, function(x){!all(is.na(x))})]
distinct_changepoint3 <- distinct_changepoint3[sapply(distinct_changepoint3, function(x){!all(is.na(x))})]

# keep only changepoints in the middle
common_changepoint4 <- common_changepoint
distinct_changepoint4 <- distinct_changepoint
p <- floor(nrow(prin_df)/2)
for(i in 1:length(common_changepoint4)){
  if(common_changepoint4[[i]]$start == 1 | common_changepoint4[[i]]$end == p){
    common_changepoint4[[i]] <- NA
  }
}
for(i in 1:length(distinct_changepoint4)){
  if(distinct_changepoint4[[i]]$start == 1 | distinct_changepoint4[[i]]$end == p){
    distinct_changepoint4[[i]] <- NA
  }
}
common_changepoint4 <- common_changepoint4[sapply(common_changepoint4, function(x){!all(is.na(x))})]
distinct_changepoint4 <- distinct_changepoint4[sapply(distinct_changepoint4, function(x){!all(is.na(x))})]

############

col_vec <- scales::hue_pal()(length(unique(new_celltype)))[c(1,4,6)]
new_celltype <- mbrain2@meta.data$celltype
wsnn_vec <- as.character(mbrain2@meta.data$wsnn_res.2)
selected_clusters <- c("17", "16", "0", "4", "2", "1", "5", "7", "12", "20")
new_celltype[!wsnn_vec %in% selected_clusters] <- NA
new_celltype[new_celltype %in% c("Forebrain GABAergic", "Glioblast", "Oligodendrocyte")] <- NA
mbrain2[["newvec"]] <- new_celltype
plot1 <- Seurat::DimPlot(mbrain2, reduction = "both", 
                         group.by = "newvec", label = TRUE,
                         repel = TRUE, label.size = 2.5, cols = col_vec)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo\nBoth Everything (Selected)") 
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both_full.png",
                plot1, device = "png", width = 5, height = 5, units = "in")


col_palette <- scales::hue_pal()(length(unique(clustering)))
col_vec <- sapply(rownames(res$adj_mat), function(x){
  col_palette[which(sort(as.numeric(unique(clustering))) == as.numeric(x))]
})

png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph_meta(mbrain2[["both"]]@cell.embeddings,
           clustering = clustering,
           adj_mat = res$adj_mat,
           col_vec = col_vec,
           cex_missing = 1,
           cex_normal = 2,
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()

png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both_diffusion_radial.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph_meta(mbrain2[["both"]]@cell.embeddings,
           clustering = clustering,
           adj_mat = res$adj_mat,
           feature_vec = dist_mat[radial_start,],
           zlim = range(dist_mat),
           cex_missing = 1,
           cex_normal = 2,
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()


png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both_diffusion_cortical.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph_meta(mbrain2[["both"]]@cell.embeddings,
           clustering = clustering,
           adj_mat = res$adj_mat,
           feature_vec = dist_mat[cortical_end,],
           zlim = range(dist_mat),
           cex_missing = 1,
           cex_normal = 2,
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()


png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_both_pseudotime.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot_graph(mbrain2[["both"]]@cell.embeddings,
           clustering = clustering,
           adj_mat = res$adj_mat,
           feature_vec = avg_rank,
           cex_missing = 1,
           cex_normal = 2,
           xlab = "both_1",
           ylab = "both_2",
           main = "SNN graph")
graphics.off()

plot1 <- plot_feature(mbrain2, 
                      cell_names = prin_df$cell,
                      feature_vec = prin_df$ord,
                      title = "Estimated pseudotime ordering")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_full_pseudotime.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

##########################

gene_vec <- c("Satb2", "Itpr1", "9130024F11Rik", "Gpm6a", "Frmd4a", "Clstn2",
              "Myt1l", "Ptprd", "Gria2", "Nav3")
gene_vec <- gene_vec[which(gene_vec %in% rownames(dcca_res$svd_1$v))]
png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_genes_cortical.png",
    height = 1500, width = 2500, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_gene_ordering(dcca_res, 
                   gene_vec = gene_vec,
                   prin_df = prin_df,
                   par_mfrow = c(3,3),
                   pch = 16,
                   xlab = "Pseudotime ordering",
                   ylab = "Expression")
graphics.off()

gene_vec <- c("8030451A03Rik", "Tmem132c", "Gm5089", "Slco1c1", "Wnt8b", 
              "Adamts19", "Ccdc80", "Plce1", "Tnc")
gene_vec <- gene_vec[which(gene_vec %in% rownames(dcca_res$svd_1$v))]
png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_genes_radial.png",
    height = 1500, width = 2500, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_gene_ordering(dcca_res, 
                   gene_vec = gene_vec,
                   prin_df = prin_df,
                   par_mfrow = c(3,3),
                   pch = 16,
                   xlab = "Pseudotime ordering",
                   ylab = "Expression")
graphics.off()


gene_vec <- c("Satb2", "Itpr1", "9130024F11Rik", "Gpm6a", "Frmd4a", "Clstn2",
              "Myt1l", "Ptprd", "Gria2", "Nav3")
gene_vec <- gene_vec[which(gene_vec %in% rownames(dcca_res$svd_1$v))]
png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_genes_cortical_common.png",
    height = 1500, width = 2500, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_gene_ordering(dcca_decomp$common_mat_1, 
                   gene_vec = gene_vec,
                   prin_df = prin_df,
                   par_mfrow = c(3,3),
                   pch = 16,
                   xlab = "Pseudotime ordering",
                   ylab = "Expression")
graphics.off()

gene_vec <- c("8030451A03Rik", "Tmem132c", "Gm5089", "Slco1c1", "Wnt8b", 
              "Adamts19", "Ccdc80", "Plce1", "Tnc")
gene_vec <- gene_vec[which(gene_vec %in% rownames(dcca_res$svd_1$v))]
png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_genes_radial_common.png",
    height = 1500, width = 2500, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_gene_ordering(dcca_decomp$common_mat_1, 
                   gene_vec = gene_vec,
                   prin_df = prin_df,
                   par_mfrow = c(3,3),
                   pch = 16,
                   xlab = "Pseudotime ordering",
                   ylab = "Expression")
graphics.off()


gene_vec <- c("Satb2", "Itpr1", "9130024F11Rik", "Gpm6a", "Frmd4a", "Clstn2",
              "Myt1l", "Ptprd", "Gria2", "Nav3")
gene_vec <- gene_vec[which(gene_vec %in% rownames(dcca_res$svd_1$v))]
png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_genes_cortical_distinct.png",
    height = 1500, width = 2500, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_gene_ordering(dcca_decomp$distinct_mat_1, 
                   gene_vec = gene_vec,
                   prin_df = prin_df,
                   par_mfrow = c(3,3),
                   pch = 16,
                   xlab = "Pseudotime ordering",
                   ylab = "Expression")
graphics.off()

gene_vec <- c("8030451A03Rik", "Tmem132c", "Gm5089", "Slco1c1", "Wnt8b", 
              "Adamts19", "Ccdc80", "Plce1", "Tnc")
gene_vec <- gene_vec[which(gene_vec %in% rownames(dcca_res$svd_1$v))]
png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_genes_radial_distinct.png",
    height = 1500, width = 2500, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_gene_ordering(dcca_decomp$distinct_mat_1, 
                   gene_vec = gene_vec,
                   prin_df = prin_df,
                   par_mfrow = c(3,3),
                   pch = 16,
                   xlab = "Pseudotime ordering",
                   ylab = "Expression")
graphics.off()

png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_de_genes.png",
    height = 1200, width = 2000, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_gene_de(common_changepoint,
             distinct_changepoint,
             p = floor(nrow(prin_df)/2),
             type = "l",
             lwd = 2,
             xlab = "Pseudotime ordering",
             ylab = "Number of DE genes",
             main1 = paste0("DE in common space\n", length(common_changepoint), 
                            " out of ", ncol(dcca_decomp$common_mat_1), " genes"),
             main2 = paste0("DE in distinct space\n", length(distinct_changepoint), 
                            " out of ", ncol(dcca_decomp$distinct_mat_1), " genes"))
graphics.off()


png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_de_genes2.png",
    height = 1200, width = 2000, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_gene_de(common_changepoint2,
             distinct_changepoint2,
             p = floor(nrow(prin_df)/2),
             type = "l",
             lwd = 2,
             xlab = "Pseudotime ordering",
             ylab = "Number of DE genes",
             main1 = paste0("DE in common space (via P-val)\n", length(common_changepoint2), 
                            " out of ", ncol(dcca_decomp$common_mat_1), " genes"),
             main2 = paste0("DE in distinct space (via P-val)\n", length(distinct_changepoint2), 
                            " out of ", ncol(dcca_decomp$distinct_mat_1), " genes"),
             cex.main = 1)
graphics.off()


png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_de_genes3.png",
    height = 1200, width = 2000, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_gene_de(common_changepoint3,
             distinct_changepoint3,
             p = floor(nrow(prin_df)/2),
             type = "l",
             lwd = 2,
             xlab = "Pseudotime ordering",
             ylab = "Number of DE genes",
             main1 = paste0("DE in common space (via length)\n", length(common_changepoint3), 
                            " out of ", ncol(dcca_decomp$common_mat_1), " genes"),
             main2 = paste0("DE in distinct space (via length)\n", length(distinct_changepoint3), 
                            " out of ", ncol(dcca_decomp$distinct_mat_1), " genes"),
             cex.main = 1)
graphics.off()


png("../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_de_genes4.png",
    height = 1200, width = 2000, units = "px", res = 300)
par(mar = c(4,4,4,0.5))
plot_gene_de(common_changepoint4,
             distinct_changepoint4,
             p = floor(nrow(prin_df)/2),
             type = "l",
             lwd = 2,
             xlab = "Pseudotime ordering",
             ylab = "Number of DE genes",
             main1 = paste0("DE in common space (middle)\n", length(common_changepoint4), 
                            " out of ", ncol(dcca_decomp$common_mat_1), " genes"),
             main2 = paste0("DE in distinct space (middle)\n", length(distinct_changepoint4), 
                            " out of ", ncol(dcca_decomp$distinct_mat_1), " genes"),
             cex.main = 1)
graphics.off()

########################

p <- floor(nrow(prin_df)/2)
intersection_genes <- intersect(names(common_changepoint), names(distinct_changepoint))
ex1 <- setdiff(names(common_changepoint), intersection_genes)
ex1_gene <- names(which.max(sapply(ex1, function(gene){
  idx <- which(names(common_changepoint) == gene)
  tmp <- common_changepoint[[idx]]$end - common_changepoint[[idx]]$start
  if(tmp < p/2 & tmp > p/10) return(common_changepoint[[idx]]$diff) else return(0)
})))

ex2 <- setdiff(names(distinct_changepoint), intersection_genes)
ex2_gene <- names(which.max(sapply(ex2, function(gene){
  idx <- which(names(distinct_changepoint) == gene)
  tmp <- distinct_changepoint[[idx]]$end - distinct_changepoint[[idx]]$start
  if(tmp < p/2 & tmp > p/10) return(distinct_changepoint[[idx]]$diff) else return(0)
})))

len_spacing <- sapply(intersection_genes, function(gene){
  idx1 <- which(names(common_changepoint) == gene)
  idx2 <- which(names(distinct_changepoint) == gene)
  
  vec1 <- common_changepoint[[idx1]]$start:common_changepoint[[idx1]]$end
  vec2 <- distinct_changepoint[[idx2]]$start:distinct_changepoint[[idx2]]$end
  
  len1 <- common_changepoint[[idx1]]$end - common_changepoint[[idx1]]$start
  len2 <- distinct_changepoint[[idx2]]$end - distinct_changepoint[[idx2]]$start
  
  if(length(intersect(vec1, vec2)) > 0 | len1 > p/2 | len1 < p/10 | len2 > p/2 | len2 < p/10){
    return(0)
  } else {
    return(min(common_changepoint[[idx1]]$diff, distinct_changepoint[[idx2]]$diff))
  }
})
ex3_gene <- intersection_genes[which.max(len_spacing)]

plot_gene_de_decomp(common_mat = dcca_decomp$common_mat_1,
                    distinct_mat = dcca_decomp$distinct_mat_1,
                    gene_vec = c("Itpr1", ex3_gene, ex1_gene),
                    prin_df = prin_df,
                    common_changepoint = common_changepoint,
                    distinct_changepoint = distinct_changepoint,
                    filepath = "../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_snn_de", 
                    height = 1200,
                    width = 1200)

