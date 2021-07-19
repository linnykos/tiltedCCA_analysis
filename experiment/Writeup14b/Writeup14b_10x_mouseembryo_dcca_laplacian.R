rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
source("../Writeup14/Writeup14_peakcalling_function.R")
date_of_run <- Sys.time(); session_info <- sessionInfo()

mat_1 <- t(mbrain[["SCT"]]@scale.data)
Seurat::DefaultAssay(mbrain) <- "ATAC"
mat_2 <- Matrix::t(mbrain[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain),])

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
metadata <- mbrain@meta.data

#######

# set.seed(10)
# mbrain <- FindMultiModalNeighbors(mbrain, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
# mbrain <- FindClusters(mbrain, graph.name = "wsnn", algorithm = 3, resolution=22,verbose = FALSE)
# ext.upstream=500000
# ext.downstream=500000
# # include some cell features object, to aggregate over.
# include.cell.features=c("nCount_RNA", "nCount_ATAC", "wsnn_res.22")
# gene_vec <- c("Pax6", "Neurog2", "Eomes", "Neurod1", "Neurod2", "Neurod4", "Neurod6", "Tbr1",
#               "Sox5", "Fezf2", "Ikzf1", "Foxg1", "Tle4", "Bcl11b", "Nr2f1", "Tbr1",
#               "Satb2", "Pou3f2", "Pou3f3")
# gene_vec <- gene_vec[gene_vec %in% colnames(mat_1)]
# myobj <- getPeaksForGenes(mbrain, gene.names = gene_vec, 
#                           include.cell.features=include.cell.features,
#                           ext.upstream=ext.upstream, ext.downstream=ext.downstream)
# myobj <- aggregateCells(myobj, aggregate.over="wsnn_res.22")
# myobj <- findLinks(myobj)

#######

cell_idx <- which(metadata$label_Savercat %in% c("Oligodendrocyte", "Hindbrain glycinergic", "Midbrain glutamatergic",
                                                 "Forebrain glutamatergic", "Radial glia", "Forebrain GABAergic", 
                                                 "Neuroblast", "Cajal-Retzius", "Mixed region GABAergic", 
                                                 "Glioblast", "Cortical or hippocampal glutamatergic"))
set.seed(10)
rank_1 <- 30; rank_2 <- 50; nn <- 15
mat_1 <- mat_1[cell_idx,]; mat_2 <- mat_2[cell_idx,]
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = T,
                                      meta_clustering = NA, num_neigh = nn, 
                                      apply_shrinkage = F, fix_distinct_perc = F, 
                                      verbose = T) 

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)
mat_1_denoised <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
mat_2_denoised <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2

############

mbrain2 <- Seurat::CreateSeuratObject(counts = t(mat_1_denoised))
mbrain2[["label_Savercat"]] <- metadata$label_Savercat[cell_idx]
membership_vec <- as.factor(metadata$label_Savercat[cell_idx])

title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                         data_1 = T, data_2 = F,
                                         radius_quantile = 0.5,
                                         bool_matrix = T, verbose = T)

#compute all rna basis vectors
k_max <- 200
c_eig <- multiomicCCA::compute_laplacian(rna_frnn$c_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                         colname_vec = paste0("clap_", 1:k_max))
d_eig <- multiomicCCA::compute_laplacian(rna_frnn$d_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                         colname_vec = paste0("dlap_", 1:k_max))
e_eig <- multiomicCCA::compute_laplacian(rna_frnn$e_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                         colname_vec = paste0("elap_", 1:k_max))

mbrain2[["clap"]] <- Seurat::CreateDimReducObject(embedding = c_eig, key = "clap", assay = "RNA")
mbrain2[["dlap"]] <- Seurat::CreateDimReducObject(embedding = d_eig, key = "dlap", assay = "RNA")
mbrain2[["elap"]] <- Seurat::CreateDimReducObject(embedding = e_eig, key = "elap", assay = "RNA")

set.seed(10)
rna_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = T, data_2 = F,
                                                 c_g = rna_frnn$c_g, d_g = rna_frnn$d_g, 
                                                 only_embedding = T, verbose = T, 
                                                 sampling_type = "adaptive_gaussian",
                                                 keep_nn = T)
mbrain2[["common"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[1]], key = "UMAP", assay = "RNA")
mbrain2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[2]], key = "UMAP", assay = "RNA")
mbrain2[["everything"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[3]], key = "UMAP", assay = "RNA")

################

set.seed(10)
atac_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                          data_1 = F, data_2 = T,
                                          radius_quantile = 0.5,
                                          bool_matrix = T, verbose = T)

#compute all the degree vectors
k_max <- 200
c_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$c_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                          colname_vec = paste0("clap_", 1:k_max))
d_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$d_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                          colname_vec = paste0("dlap_", 1:k_max))
e_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$e_g, k_max = k_max, rowname_vec = colnames(mbrain2), 
                                          colname_vec = paste0("elap_", 1:k_max))

mbrain2[["clap2"]] <- Seurat::CreateDimReducObject(embedding = c_eig2, key = "clap", assay = "RNA")
mbrain2[["dlap2"]] <- Seurat::CreateDimReducObject(embedding = d_eig2, key = "dlap", assay = "RNA")
mbrain2[["elap2"]] <- Seurat::CreateDimReducObject(embedding = e_eig2, key = "elap", assay = "RNA")

set.seed(10)
atac_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = F, data_2 = T,
                                                  c_g = atac_frnn$c_g, d_g = atac_frnn$d_g, 
                                                  only_embedding = T, verbose = T, 
                                                  sampling_type = "adaptive_gaussian",
                                                  keep_nn = T)
mbrain2[["common2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[1]], key = "UMAP", assay = "RNA")
mbrain2[["distinct2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[2]], key = "UMAP", assay = "RNA")
mbrain2[["everything2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[3]], key = "UMAP", assay = "RNA")

##########################

set.seed(10)
combined_g <- multiomicCCA::combine_frnn(dcca_res, g_1 = rna_frnn$c_g,
                                          g_2 = atac_frnn$c_g, nn = 15)
combined_g <- SeuratObject::as.Graph(combined_g)
set.seed(10)
combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings
mbrain2[["combined"]] <- Seurat::CreateDimReducObject(embedding = combined_common_umap, key = "UMAP", assay = "RNA")

set.seed(10)
both_embeddings <- multiomicCCA::plot_embeddings(dcca_res, membership_vec = membership_vec,
                                                 data_1 = T, data_2 = T, add_noise = F,
                                                 only_embedding = T)
mbrain2[["both"]] <- Seurat::CreateDimReducObject(embedding = both_embeddings$everything, key = "UMAP", assay = "RNA")

##########################################3
# done with all the calculations (for now). now plot

# plot RNA embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(mbrain2, reduction = main_vec[i], 
                           group.by = "label_Savercat", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (RNA)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_rna_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot ATAC embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(mbrain2, reduction = paste0(main_vec[i], "2"), 
                           group.by = "label_Savercat", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo (ATAC)\n", title_vec[i])) 
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_atac_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot the combined view
plot1 <- Seurat::DimPlot(mbrain2, reduction = "combined", 
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo\nBoth Common") 
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_both_common_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

# plot the everything combined for both
plot1 <- Seurat::DimPlot(mbrain2, reduction = "both", 
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo\nBoth Everything") 
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_both_everything_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


#######################################
# now plot the bases -- first RNA
plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("clap_", 1:16), reduction = "combined")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_rna_common_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("dlap_", 1:16), reduction = "distinct")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_rna_distinct_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("elap_", 1:16), reduction = "both")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_rna_everything_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

# next ATAC
plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("clap2_", 1:16), reduction = "combined")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_atac_common_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("dlap2_", 1:16), reduction = "distinct2")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_atac_distinct_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

plot1 <- Seurat::FeaturePlot(mbrain2, features = paste0("elap2_", 1:16), reduction = "both")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_atac_everything_basis.png"),
                plot1, device = "png", width = 16, height = 12, units = "in")

##################

# compute local enrichment
set.seed(10)
rna_local <- multiomicCCA::clisi_information(rna_frnn$c_g, rna_frnn$d_g, rna_frnn$e_g, 
                                             membership_vec = membership_vec)
rna_local$common_clisi$membership_info
rna_local$distinct_clisi$membership_info

set.seed(10)
atac_local <- multiomicCCA::clisi_information(atac_frnn$c_g, atac_frnn$d_g, atac_frnn$e_g, 
                                             membership_vec = membership_vec)
atac_local$common_clisi$membership_info
atac_local$distinct_clisi$membership_info

n <- nrow(rna_local$common_clisi$cell_info)
k <- nrow(rna_local$common_clisi$membership_info)
df <- data.frame(celltype = as.factor(c(paste0(as.character(rna_local$common_clisi$cell_info$celltype), "0"), 
                                        as.character(rna_local$common_clisi$membership_info$celltype))), 
                 common = c(rna_local$common_clisi$cell_info$clisi_score, rna_local$common_clisi$membership_info$mean_clisi),
                 distinct = c(rna_local$distinct_clisi$cell_info$clisi_score, rna_local$distinct_clisi$membership_info$mean_clisi),
                 category = as.factor(c(rep(0, n), rep(1, k))))
col_vec <- scales::hue_pal()(k)
bg_col_vec <- multiomicCCA:::.adjust_colors(col_vec, l_bg = 75, c_bg = 50, alpha_bg = 0.5)
all_col_vec <- c(col_vec, bg_col_vec)
tmp <- rna_local$common_clisi$membership_info$celltype
names(all_col_vec) <- c(tmp, paste0(tmp, "0"))
custom_colors <- ggplot2::scale_colour_manual(values = all_col_vec)

plot1 <- ggplot2::ggplot(data = subset(df, category == 0), ggplot2::aes(x = distinct, y = common, color = celltype))
plot1 <- plot1 + ggplot2::geom_point()
plot1 <- plot1 + ggplot2::xlim(1, 0) + ggplot2::ylim(0, 1)
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, category == 1), 
                                     ggplot2::aes(x = distinct, y = common), 
                                     size = 3, color = "black")
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, category == 1), 
                                     ggplot2::aes(x = distinct, y = common), 
                                     size = 2.5, color = "white")
plot1 <- plot1 + ggplot2::geom_point(data = subset(df, category == 1), 
                                     ggplot2::aes(x = distinct, y = common,
                                                  color = celltype), 
                                     size = 2)
plot1 <- plot1 + custom_colors
plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, category == 1), ggplot2::aes(label = celltype),
                                          color = "black",
                                          segment.color = "grey50",
                                          size = 2)
# see directions at https://ggrepel.slowkow.com/articles/examples.html
plot1 <- plot1 + ggplot2::xlab("Distinct enrichment")
plot1 <- plot1 + ggplot2::ylab("Common enrichment")
plot1 <- plot1 + ggplot2::ggtitle("RNA")
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_rna_enrichment.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")


plot1 <- plot1 + ggplot2::scale_size_continuous(range = c(1, 2))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c("black", "green", "red"))
plot1 <- plot1 + ggplot2::geom_hline(yintercept = 0, linetype="dashed",
               color = "orange")




########################
########################
########################

# p1 <- ncol(mat_1); p2 <- ncol(mat_2)
# gene_smoothed <- lapply(1:p1, function(j){
#   if(j %% floor(p1/10) == 0) cat('*')
#   
#   c_res <- compute_smooth_signal(mat_1_denoised[,j], c_eig)
#   d_res <- compute_smooth_signal(mat_1_denoised[,j], d_eig)
#   e_res <- compute_smooth_signal(mat_1_denoised[,j], e_eig)
#   
#   list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
#        d_variance = d_res$variance, d_r2 = d_res$r_squared,
#        e_variance = e_res$variance, e_r2 = e_res$r_squared)
# })
# 
# # plot some of the genes with smallest r2
# idx <- order(sapply(gene_smoothed, function(x){max(x$c_r2, x$d_r2)}), decreasing = F)[1:5]
# for(j in 1:length(idx)){
#   c_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], c_eig)
#   d_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], d_eig)
#   e_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], e_eig)
#   
#   multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_1)[idx[j]], 
#                                prefix = "RNA", e_vec = mat_1_denoised[,idx[j]],
#                                c_vec = dcca_decomp$common_mat_1[,idx[j]],
#                                d_vec = dcca_decomp$distinct_mat_1[,idx[j]],
#                                e_res = e_res, c_res = c_res, d_res = d_res,
#                                filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_lowesetr2_gene", j, ".png"))
# }
# 
# # plot some of the genes with highest r2
# idx <- order(sapply(gene_smoothed, function(x){max(x$c_r2, x$d_r2)}), decreasing = T)[1:5]
# for(j in 1:length(idx)){
#   c_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], c_eig)
#   d_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], d_eig)
#   e_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], e_eig)
#   
#   multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_1)[idx[j]], 
#                                prefix = "RNA", e_vec = mat_1_denoised[,idx[j]],
#                                c_vec = dcca_decomp$common_mat_1[,idx[j]],
#                                d_vec = dcca_decomp$distinct_mat_1[,idx[j]],
#                                e_res = e_res, c_res = c_res, d_res = d_res,
#                                filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_highestr2_gene", j, ".png"))
# }

# # plot some of the genes with biggest difference
# tmp <- sapply(gene_smoothed, function(x){
#   (x$d_variance - x$c_variance)/x$e_variance
# })
# idx <- order(tmp, decreasing = T)[1:5]
# for(j in 1:length(idx)){
#   c_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], c_eig)
#   d_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], d_eig)
#   e_res <- compute_smooth_signal(mat_1_denoised[,idx[j]], e_eig)
#   
#   multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_1)[idx[j]], 
#                                prefix = "RNA", e_vec = mat_1_denoised[,idx[j]],
#                                c_vec = dcca_decomp$common_mat_1[,idx[j]],
#                                d_vec = dcca_decomp$distinct_mat_1[,idx[j]],
#                                e_res = e_res, c_res = c_res, d_res = d_res,
#                                filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_biggestdiff_gene", j, ".png"))
# }

# # make a ggplot of the distinct-commmon
# vec1 <- sapply(gene_smoothed, function(x){
#   (x$d_variance - x$c_variance)/x$e_variance
# })
# vec2 <- sapply(gene_smoothed, function(x){x$d_variance - x$c_variance})
# gene_factor <- rep(0, p1); gene_factor[order(vec2, decreasing = T)[1:100]] <- 1
# gene_factor[which(colnames(mat_1_denoised) %in% gene_vec)] <- 2
# size_vec <- rep(1, p1); size_vec[which(colnames(mat_1_denoised) %in% gene_vec)] <- 2
# df <- data.frame(val1 = vec1, val2 = vec2,
#                  idx = rank(vec1), 
#                  gene_factor = gene_factor,
#                  size_vec = size_vec,
#                  text = colnames(mat_1_denoised))
# # reshuffle
# df <- df[order(df$gene_factor, decreasing = F),]
# df$gene_factor <- as.factor(df$gene_factor)
# plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = idx, y = val1)) 
# plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = gene_factor, size = size_vec))
# plot1 <- plot1 + ggplot2::scale_size_continuous(range = c(1, 2))
# plot1 <- plot1 + ggplot2::scale_colour_manual(values=c("black", "green", "red"))
# plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, gene_factor == 2), ggplot2::aes(label = text, color = gene_factor),
#                                     box.padding = ggplot2::unit(0.5, 'lines'),
#                                     point.padding = ggplot2::unit(1.6, 'lines'),
#                                     nudge_y = 0.25, size = 2)
# plot1 <- plot1 + ggplot2::geom_hline(yintercept = 0, linetype="dashed", 
#                 color = "orange")
# plot1 <- plot1 + ggplot2::xlab("Order of genes")
# plot1 <- plot1 + ggplot2::ylab("Distinct variance - common variance")
# plot1 <- plot1 + ggplot2::ggtitle("Gene-level decomposition")
# plot1 <- plot1 + Seurat::NoLegend()
# ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_rna_gene_decomposition.png"),
#                 plot1, device = "png", width = 4, height = 4, units = "in")

# # make a ggplot of the r2
# vec1 <- sapply(gene_smoothed, function(x){min(abs(x$c_r2), abs(x$d_r2))})
# vec2 <- sapply(gene_smoothed, function(x){x$d_variance - x$c_variance})
# gene_factor <- rep(0, p1); gene_factor[order(vec2, decreasing = T)[1:100]] <- 1
# gene_factor[which(colnames(mat_1_denoised) %in% gene_vec)] <- 2
# size_vec <- rep(1, p1); size_vec[which(colnames(mat_1_denoised) %in% gene_vec)] <- 2
# df <- data.frame(val1 = vec1, val2 = vec2,
#                  idx = rank(vec1), 
#                  gene_factor = gene_factor,
#                  size_vec = size_vec,
#                  text = colnames(mat_1_denoised))
# # reshuffle
# df <- df[order(df$gene_factor, decreasing = F),]
# df$gene_factor <- as.factor(df$gene_factor)
# plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = idx, y = val1)) 
# plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = gene_factor, size = size_vec))
# plot1 <- plot1 + ggplot2::scale_size_continuous(range = c(1, 2))
# plot1 <- plot1 + ggplot2::scale_colour_manual(values=c("black", "green", "red"))
# plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, gene_factor == 2), ggplot2::aes(label = text, color = gene_factor),
#                                           box.padding = ggplot2::unit(0.5, 'lines'),
#                                           point.padding = ggplot2::unit(1.6, 'lines'),
#                                           nudge_y = -0.25, size = 2)
# plot1 <- plot1 + ggplot2::xlab("Order of genes")
# plot1 <- plot1 + ggplot2::ylab("Minimum R2 of conformity to graph structure")
# plot1 <- plot1 + ggplot2::ggtitle("Gene-level decomposition")
# plot1 <- plot1 + Seurat::NoLegend()
# ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_rna_gene_conformity.png"),
#                 plot1, device = "png", width = 4, height = 4, units = "in")


##########################################
##########################################
##########################################




#################################3
# # segway to test whether or not our way was good
# embedding <- multiomicCCA:::.prepare_embeddings(dcca_res, data_1 = F, data_2 = T, 
#                                  add_noise = F)
# set.seed(10)
# zz <- Seurat::RunUMAP(embedding[[1]])
# mbrain2[["common"]] <- Seurat::CreateDimReducObject(embedding = zz@cell.embeddings, key = "UMAP", assay = "RNA")
# set.seed(10)
# zz <- Seurat::RunUMAP(embedding[[2]])
# mbrain2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = zz@cell.embeddings, key = "UMAP", assay = "RNA")
# 
# for(i in 1:3){
#   plot1 <- Seurat::DimPlot(mbrain2, reduction = main_vec[i], 
#                            group.by = "label_Savercat", label = TRUE,
#                            repel = TRUE, label.size = 2.5)
#   plot1 <- plot1 + ggplot2::ggtitle("Mouse Embryo (ATAC)\n", title_vec[i]) 
#   plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
#   ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_dcca_atac_", main_vec[i], "_umap.png"),
#                   plot1, device = "png", width = 6, height = 5, units = "in")
# }

#################################3

####################################

# atac_smoothed <- lapply(1:p2, function(j){
#   print(j)
#   
#   c_res <- compute_smooth_signal(mat_2_denoised[,j], c_eig2)
#   d_res <- compute_smooth_signal(mat_2_denoised[,j], d_eig2)
#   e_res <- compute_smooth_signal(mat_2_denoised[,j], e_eig2)
#   
#   list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
#        d_variance = d_res$variance, d_r2 = d_res$r_squared,
#        e_variance = e_res$variance, e_r2 = e_res$r_squared)
# })
# 
# save(atac_smoothed, file = "../../../../out/Writeup14b/atac_smoothed.RData")
# load("../../../../out/Writeup14b/atac_smoothed.RData")
# 
# set.seed(10)
# atac_embeddings <- multiomicCCA::plot_embeddings(dcca_res, membership_vec, data_1 = F, data_2 = T, 
#                                                 add_noise = F, pca = F, only_embedding = T, verbose = T)
# mbrain2 <- Seurat::CreateSeuratObject(counts = t(mat_1_denoised))
# mbrain2[["label_Savercat"]] <- metadata$label_Savercat[cell_idx]
# mbrain2[["common"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[1]], key = "UMAP", assay = "RNA")
# mbrain2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[2]], key = "UMAP", assay = "RNA")
# mbrain2[["everything"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[3]], key = "UMAP", assay = "RNA")
# 
# 
# # plot some of the genes with smallest r2
# idx <- order(sapply(atac_smoothed, function(x){max(x$c_r2, x$d_r2)}), decreasing = F)[1:5]
# for(j in 1:length(idx)){
#   c_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], c_eig)
#   d_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], d_eig)
#   e_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], e_eig)
#   
#   multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_2)[idx[j]], 
#                                prefix = "ATAC", e_vec = mat_2_denoised[,idx[j]],
#                                c_vec = dcca_decomp$common_mat_2[,idx[j]],
#                                d_vec = dcca_decomp$distinct_mat_2[,idx[j]],
#                                e_res = e_res, c_res = c_res, d_res = d_res,
#                                filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_lowesetr2_atac", j, ".png"))
# }
# 
# # plot some of the genes with highest r2
# idx <- order(sapply(atac_smoothed, function(x){max(x$c_r2, x$d_r2)}), decreasing = T)[1:5]
# for(j in 1:length(idx)){
#   c_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], c_eig)
#   d_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], d_eig)
#   e_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], e_eig)
#   
#   multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_2)[idx[j]], 
#                                prefix = "ATAC", e_vec = mat_2_denoised[,idx[j]],
#                                c_vec = dcca_decomp$common_mat_2[,idx[j]],
#                                d_vec = dcca_decomp$distinct_mat_2[,idx[j]],
#                                e_res = e_res, c_res = c_res, d_res = d_res,
#                                filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_highestr2_atac", j, ".png"))
# }
# 
# # plot some of the genes with biggest difference
# tmp <- sapply(atac_smoothed, function(x){
#   max_val <- max(x$c_variance,x$d_variance)
#   min_val <- min(x$c_variance,x$d_variance)
#   (max_val - min_val)/min_val
# })
# idx <- order(tmp, decreasing = T)[1:5]
# for(j in 1:length(idx)){
#   c_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], c_eig)
#   d_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], d_eig)
#   e_res <- compute_smooth_signal(mat_2_denoised[,idx[j]], e_eig)
#   
#   multiomicCCA::plot_laplacian(mbrain2, var_name = colnames(mat_2)[idx[j]], 
#                                prefix = "ATAC", e_vec = mat_2_denoised[,idx[j]],
#                                c_vec = dcca_decomp$common_mat_2[,idx[j]],
#                                d_vec = dcca_decomp$distinct_mat_2[,idx[j]],
#                                e_res = e_res, c_res = c_res, d_res = d_res,
#                                filename = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_biggestdiff_atac", j, ".png"))
# }
# 
# 
# ########################
# 
# png(file = paste0("../../../../out/figures/Writeup14b/Writeup14b_10x_mouseembryo_laplacian_r2.png"),
#     height = 1500, width = 3000, res = 300, units = "px")
# par(mfrow = c(1,2))
# tmp <- sapply(gene_smoothed, function(x){max(x$c_r2, x$d_r2)})
# idx <- which(colnames(mat_1) %in% gene_vec)
# col_vec <- rep(rgb(0.5,0.5,0.5,0.2), ncol(mat_1)); col_vec[idx] <- "red"
# order_idx <- order(tmp, decreasing = F)
# plot(tmp[order_idx], pch = 16, col = col_vec[order_idx], xlab = "Order of genes", ylab = "R2 b/w denoised vector and frNN Laplacian",
#      main = "RNA")
# 
# idx <- unlist(lapply(1:length(gene_vec), function(i){
#   atac_names <- colnames(myobj$X.aggr[[i]]$counts)
#   include_bool <- which(myobj$X.aggr[[i]]$peaks.gr$pval.spearman < 0.05)
#   atac_names <- atac_names[include_bool]
#   atac_idx <- which(colnames(mat_2) %in% atac_names)
#   atac_idx
# }))
# tmp <- sapply(atac_smoothed, function(x){max(x$c_r2, x$d_r2)})
# col_vec <- rep(rgb(0.5,0.5,0.5,0.2), ncol(mat_1)); col_vec[idx] <- "red"
# order_idx <- order(tmp, decreasing = F)
# plot(tmp[order_idx], pch = 16, col = col_vec[order_idx], xlab = "Order of peaks", ylab = "R2 b/w denoised vector and frNN Laplacian",
#      main = "ATAC")
# graphics.off()

##################



