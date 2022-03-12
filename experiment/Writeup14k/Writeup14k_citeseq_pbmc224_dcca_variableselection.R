rm(list=ls())
library(Seurat)
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca2.RData")
load("../../../../out/Writeup14i/Writeup14i_citeseq_pbmc224_dcca_variable_full_de.RData")
dcca_res <- dcca_res2
source("../Writeup14f/gene_exploration.R")
source("variable_distinct_selection.R")

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res)

summary_mat <- compute_variable_summary(mat = mat_2b, 
                                        common_mat = dcca_decomp$common_mat_2,
                                        metacell_clustering = NA)
summary_mat <- data.frame(summary_mat)
summary_mat <- cbind(summary_mat, NA)
colnames(summary_mat) <- c("r_squared", "p_value", "celltype")
uniq_celltype <- as.character(sort(unique(pbmc$celltype.l2)))
names(de_list) <- uniq_celltype

for(i in 1:nrow(summary_mat)){
  gene_name <- rownames(summary_mat)[i]
  
  max_celltype_idx <- max(combn_mat)
  celltype_pval <- sapply(1:max_celltype_idx, function(celltype){
    idx <- which(combn_mat == celltype, arr.ind = T)[,2]
    stats::quantile(sapply(idx, function(j){
      idx <- which(rownames(de_list[[j]]) == gene_name)
      if(length(idx) == 0) return(1)
      de_list[[j]][idx, "p_val"]
    }), probs = 0.75)
  })
  
  summary_mat[i,"p_value"] <- min(celltype_pval)
  summary_mat[i,"celltype"] <- uniq_celltype[which.min(celltype_pval)]
}
summary_mat[,"p_value"] <- -log10(summary_mat[,"p_value"])
summary_mat[order(rownames(summary_mat)),]

######################

png("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_adt_exploration.png",
    height = 1250, width = 1250, res = 300, units = "px")
tmp <- summary_mat[,"p_value"]; tmp <- tmp[tmp > 0]
xmin <- floor(log10(min(tmp)))
summary_mat[which(summary_mat[,"p_value"] == 0),"p_value"] <- 10^xmin
xmax <- ceiling(max(log10(summary_mat[,"p_value"])))
plot(NA,
     xlim = c(xmin, xmax),
     ylim = c(0,1),
     bty = "n",
     xaxt = "n",
     xlab = "Separability (Cell-type -Log10(p-value))",
     ylab = "Alignment w/ common space (R^2)",
     main = "CITE-Seq: PBMC (Protein)")
tmp <- seq(xmin, xmax, by = 1)
xaxt_vec <- sort(unique(unlist(sapply(1:(length(tmp)-1), function(i){
  lower <- 10^(tmp[i]); upper <- 10^(tmp[i+1])
  seq_vec <- seq(lower, upper, length.out = 10)
  log10(seq_vec)
}))))
axis(1, at = xaxt_vec, labels = rep("", length(xaxt_vec)))
axis(1, at = tmp, labels = as.character(10^tmp))
for(y_val in seq(0,1,by=0.1)){
  lines(c(-1e5,1e5), rep(y_val,2), lty = 2, col = "gray", lwd = 1)
}
tmp <- pmax(seq(0,length(xaxt_vec), by = 2),1)
major <- c(1, tmp[tmp %% 10 == 0])
minor <- tmp[!tmp %in% major]
for(x_val in xaxt_vec[minor]){
  lines(rep(x_val,2), c(-1e5,1e5), lty = 3, col = "gray", lwd = 0.5)
}
for(x_val in xaxt_vec[major]){
  lines(rep(x_val,2), c(-1e5,1e5), lty = 2, col = "gray", lwd = 1.5)
}

points(log10(summary_mat[,"p_value"]), summary_mat[,"r_squared"], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 1))
graphics.off()


########### ## ggrepel version of the plot
df <- as.data.frame(summary_mat)
df$Name <- rownames(df)
p1 <- ggplot2::ggplot(df, ggplot2::aes(x = p_value, y = r_squared))
p1 <- p1 + ggplot2::geom_point()
p1 <- p1 + ggplot2::scale_x_log10()
# p1 <- p1 + ggplot2::xlim(1, max(df$kl_div))
p1 <- p1 + ggplot2::ylim(0, 1)
p1 <- p1 + ggrepel::geom_text_repel(data = df, ggplot2::aes(label = Name))
p1 <- p1 + Seurat::NoLegend()
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_adt_exploration2.png",
                p1, device = "png", width = 5, height = 5, units = "in")

############################

antibodies_dat <- read.csv("../../../../data/biolegend_antibodies.csv", header = T)
celltype_names <- c("B Cells", 
                    "T Cells", 
                    "Platelets",
                    "Natural Killer Cells",
                    "Hematopoietic Stem Cells",
                    "Monocyte",
                    "Plasmacytoid Dendritic Cells",
                    "Dendritic Cells")
antibody_list <- sapply(celltype_names, function(celltype){
  idx <- grep(celltype, antibodies_dat[,1], fixed = T)
  antibody_vec <- antibodies_dat[idx,2]
  col_idx <- unlist(sapply(antibody_vec, function(antibody){
    idx1 <- which(colnames(mat_2) == antibody)
    idx2 <- grep(paste0(antibody, "-"), colnames(mat_2))
    unique(c(idx1, idx2))
  }))
  colnames(mat_2)[col_idx]
})
names(antibody_list) <- sapply(names(antibody_list), function(str){
  sub(" ", "-", str)
})


for(i in 1:length(antibody_list)){
  df <- as.data.frame(summary_mat)
  df$Labeling <- rep(0, nrow(df))
  df$Labeling[which(rownames(df) %in% antibody_list[[i]])] <- 1
  df$Labeling <- as.factor(df$Labeling)
  df$Name <- rownames(df)
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = p_value, y = r_squared))
  p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = Labeling))
  p1 <- p1 + ggplot2::scale_color_manual(values = c("black", "red"))
  p1 <- p1 + ggplot2::scale_x_log10()
  p1 <- p1 + ggplot2::ylim(0, 1)
  p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, Labeling == 1), 
                                      ggplot2::aes(label = Name, color = Labeling))
  p1 <- p1 + Seurat::NoLegend()
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_adt_exploration_", names(antibody_list)[i], ".png"),
                  p1, device = "png", width = 5, height = 5, units = "in")
}

#################################

var_selection_res <- variable_distinct_selection(dcca_res2,
                                                 significance_vec = summary_mat[,"p_value"], #larger is more significant
                                                 max_variables = 10,
                                                 cor_threshold = 0.9,
                                                 verbose = T)
antibody_vec <- var_selection_res$selected_variables
summary_mat[antibody_vec,]

df <- as.data.frame(summary_mat)
df$Labeling <- rep(0, nrow(df))
df$Labeling[which(rownames(df) %in% antibody_vec)] <- 1
df$Labeling <- as.factor(df$Labeling)
df$Name <- rownames(df)
p1 <- ggplot2::ggplot(df, ggplot2::aes(x = p_value, y = r_squared))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = Labeling))
p1 <- p1 + ggplot2::scale_color_manual(values = c("black", "red"))
p1 <- p1 + ggplot2::scale_x_log10()
p1 <- p1 + ggplot2::ylim(0, 1)
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, Labeling == 1), 
                                    ggplot2::aes(label = Name, color = Labeling))
p1 <- p1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_adt_exploration_selected.png"),
                p1, device = "png", width = 5, height = 5, units = "in")


anchor_name <- "rna.umap"
other_names <- c("adt.umap")

for(umap_name in other_names){
  print(umap_name)
  u_mat1 <- pbmc[[anchor_name]]@cell.embeddings
  u_mat2 <- pbmc[[umap_name]]@cell.embeddings
  tmp <- svd(t(u_mat1) %*% u_mat2)
  rotation_mat <- tmp$u %*% t(tmp$v)
  tmp <- u_mat2 %*% t(rotation_mat)
  rownames(tmp) <- rownames(pbmc@meta.data)
  colnames(tmp) <- colnames(pbmc[[umap_name]]@cell.embeddings)
  pbmc[[umap_name]]@cell.embeddings <- tmp
}

plot1 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\n(ADT)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_dcca_adtUmap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

nrow_val <- ceiling(sqrt(length(antibody_vec)))
ncol_val <- ceiling(length(antibody_vec)/nrow_val)
plot1 <- Seurat::FeaturePlot(pbmc, 
                             reduction = "adt.umap",
                             slot = "scale.data",
                             features = antibody_vec,
                             ncol = ncol_val)

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_dcca_adtUmap_selected.png"),
                plot1, device = "png", width = 4*ncol_val, height = 4*nrow_val, units = "in")

##########################################

set.seed(10)
pbmc2 <- pbmc
pbmc2[["adt.umap"]] <- NULL
pbmc2[["apca"]] <- NULL

# compute new ADT PCA
print("Computing ADT PCA")
Seurat::DefaultAssay(pbmc2) <- "ADT"
pbmc2[["ADT"]]@scale.data <- pbmc2[["ADT"]]@scale.data[antibody_vec,]

print(dim(pbmc2[["ADT"]]@scale.data))
pbmc2 <- Seurat::RunPCA(
  pbmc2, 
  reduction.name = 'apca',
  features = rownames(pbmc2[["ADT"]]@scale.data),
  verbose = F
)

# compute new ADT umap
print("Computing ADT UMAP")
set.seed(10)
pbmc2 <- Seurat::RunUMAP(
  pbmc2, 
  reduction = 'apca', 
  dims = 1:9, assay = 'ADT', 
  reduction.name = 'adt.umap', 
  reduction.key = 'adtUMAP_'
)

# compute new multimodal neighbor
print("Computing WNN")
pbmc2 <- Seurat::FindMultiModalNeighbors(
  pbmc2, 
  reduction.list = list("pca", "apca"), 
  dims.list = list(1:rank_1, 1:9), 
  modality.weight.name = "RNA.weight"
)

# compute new WNN UMAP
print("Computing WNN UMAP")
set.seed(10)
pbmc2 <- Seurat::RunUMAP(
  pbmc2, 
  nn.name = "weighted.nn", 
  reduction.name = "wnn.umap", 
  reduction.key = "wnnUMAP_"
)

# compute consensus PCA's UMAP
print("Computing Consensus PCA")
n <- ncol(pbmc2)
svd_1 <- svd(pbmc2[["pca"]]@cell.embeddings)
embedding_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
svd_2 <- svd(pbmc2[["apca"]]@cell.embeddings)
embedding_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))
embedding_all <- cbind(embedding_1, embedding_2)
rownames(embedding_all) <- colnames(pbmc2)
pca_res <- stats::prcomp(embedding_all, center = TRUE, scale. = FALSE)
consensus_mat <- pca_res$x[,1:min(rank_1, rank_2)] 

set.seed(10)
consensus_umap <- Seurat::RunUMAP(consensus_mat, 
                                  reduction.key = "umapConsensusPCA_")
pbmc2[["consensus.umap"]] <- Seurat::CreateDimReducObject(consensus_umap@cell.embeddings,
                                                          assay = "ADT")

# rotate all the UMAPs
print("Rotating UMAPs")
anchor_name <- "rna.umap"
other_names <- c("adt.umap", "wnn.umap", "consensus.umap")

for(umap_name in other_names){
  print(umap_name)
  u_mat1 <- pbmc2[[anchor_name]]@cell.embeddings
  u_mat2 <- pbmc2[[umap_name]]@cell.embeddings
  tmp <- svd(t(u_mat1) %*% u_mat2)
  rotation_mat <- tmp$u %*% t(tmp$v)
  tmp <- u_mat2 %*% t(rotation_mat)
  rownames(tmp) <- rownames(pbmc2@meta.data)
  colnames(tmp) <- colnames(pbmc2[[umap_name]]@cell.embeddings)
  pbmc2[[umap_name]]@cell.embeddings <- tmp
}

# make all the plots
print("Plotting")
reduction_vec <- c(other_names)
group_vec <- c("celltype.l2")
main_vec <- c("(ADT)", "(WNN)", "(Consensus PCA)")
file_vec <- c("adt", "wnn", "consensuspca")

for(i in 1:length(reduction_vec)){
  for(j in 1:length(group_vec)){
    plot1 <- Seurat::DimPlot(pbmc2, reduction = reduction_vec[i],
                             group.by = group_vec[j], label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\n", main_vec[i]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_dcca_selected_", file_vec[i], ".png"),
                    plot1, device = "png", width = 6, height = 5, units = "in")
  }
}

##################################

adt_mat <- tcrossprod(multiomicCCA:::.mult_mat_vec(dcca_res2$svd_2$u, dcca_res2$svd_2$d), dcca_res2$svd_2$v)
cor_mat <- stats::cor(adt_mat)
cor_mat2 <- abs(cor_mat)^4
diag(cor_mat2) <- 0
library(igraph)
# https://stackoverflow.com/questions/21300821/igraph-creating-a-weighted-adjacency-matrix
g <- igraph::graph.adjacency(cor_mat2, mode="undirected", weighted=TRUE)
igraph::E(g)$width <- igraph::E(g)$weight

png("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_adt_correlation_graph.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = rep(0.5,4))
set.seed(10)
plot(g, vertex.size = 3, vertex.label = "")
graphics.off()

x_vec <- summary_mat[,"r_squared"]^4
color_palette <- grDevices::colorRampPalette(c("white", 2))(20)
seq_vec <- seq(min(x_vec), max(x_vec), length.out = 20)
col_vec <- sapply(names(igraph::V(g)), function(i){
  val <- x_vec[i]
  idx <- which.min(abs(seq_vec - val))
  color_palette[idx]
})
igraph::V(g)$color <- col_vec
png("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_adt_correlation_graph_rsquared.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = rep(0.5,4))
set.seed(10)
plot(g, vertex.size = 3, vertex.label = "", vertex.color = igraph::V(g)$color )
graphics.off()

x_vec <- summary_mat[,"p_value"]
color_palette <- grDevices::colorRampPalette(c("white", 2))(20)
seq_vec <- seq(min(x_vec), max(x_vec), length.out = 20)
col_vec <- sapply(names(igraph::V(g)), function(i){
  val <- x_vec[i]
  idx <- which.min(abs(seq_vec - val))
  color_palette[idx]
})
igraph::V(g)$color <- col_vec
png("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_adt_correlation_graph_pvalue.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = rep(0.5,4))
set.seed(10)
plot(g, vertex.size = 3, vertex.label = "", vertex.color = igraph::V(g)$color )
graphics.off()

col_vec <- sapply(names(igraph::V(g)), function(i){
  if(i %in% antibody_vec) "deepskyblue2" else "white"
})
igraph::V(g)$color <- col_vec
png("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_adt_correlation_graph_selected.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = rep(0.5,4))
set.seed(10)
plot(g, vertex.size = 3, vertex.label = "", vertex.color = igraph::V(g)$color )
graphics.off()

####

x_vec <- summary_mat[,"r_squared"]^4
color_palette <- grDevices::colorRampPalette(c("white", 2))(20)
seq_vec <- seq(min(x_vec), max(x_vec), length.out = 20)
col_vec_r2 <- sapply(names(igraph::V(g)), function(i){
  val <- x_vec[i]
  idx <- which.min(abs(seq_vec - val))
  color_palette[idx]
})
x_vec <- summary_mat[,"p_value"]
color_palette <- grDevices::colorRampPalette(c("white", 2))(20)
seq_vec <- seq(min(x_vec), max(x_vec), length.out = 20)
col_vec_pval <- sapply(names(igraph::V(g)), function(i){
  val <- x_vec[i]
  idx <- which.min(abs(seq_vec - val))
  color_palette[idx]
})

for(iter in 1:4){
  col_vec_r2_tmp <- col_vec_r2; col_vec_pval_tmp <- col_vec_pval
  col_vec_r2_tmp[which(!names(col_vec_r2) %in% var_selection_res$candidate_list[[iter]])]<- "gray"
  col_vec_pval_tmp[which(!names(col_vec_pval) %in% var_selection_res$candidate_list[[iter]])]<- "gray"
  
  igraph::V(g)$color <- col_vec_r2_tmp
  png(paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_adt_correlation_graph_rsquared_iter", iter, ".png"),
      height = 2500, width = 2500, units = "px", res = 300)
  par(mar = rep(0.5,4))
  set.seed(10)
  plot(g, vertex.size = 3, vertex.label = "", vertex.color = igraph::V(g)$color )
  graphics.off()
  
  igraph::V(g)$color <- col_vec_pval_tmp
  png(paste0("../../../../out/figures/Writeup14k/Writeup14k_citeseq_pbmc224_adt_correlation_graph_pvalue_iter", iter, ".png"),
      height = 2500, width = 2500, units = "px", res = 300)
  par(mar = rep(0.5,4))
  set.seed(10)
  plot(g, vertex.size = 3, vertex.label = "", vertex.color = igraph::V(g)$color )
  graphics.off()
}