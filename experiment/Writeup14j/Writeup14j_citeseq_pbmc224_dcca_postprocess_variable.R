rm(list=ls())
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca.RData")
load("../../../../out/Writeup14i/Writeup14i_citeseq_pbmc224_dcca_variable_full_de.RData")
source("../Writeup14f/gene_exploration.R")
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca_tmp.RData")
dcca_res <- res
class(dcca_res) <- "dcca"
library(Seurat)

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res)

summary_mat <- compute_variable_summary(mat = mat_2b, 
                                        common_mat = dcca_decomp$common_mat_2,
                                        metacell_clustering = NA)
for(i in 1:nrow(summary_mat)){
  gene_name <- rownames(summary_mat)[i]
  
  val <- stats::median(sapply(1:length(de_list), function(j){
    idx <- which(rownames(de_list[[j]]) == gene_name)
    if(length(idx) == 0) return(1)
    de_list[[j]][idx, "p_val"]
  }))
  
  summary_mat[i,"kl_div"] <- val
}
summary_mat[,"kl_div"] <- -log10(summary_mat[,"kl_div"])

png("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_adt_exploration.png",
    height = 1250, width = 1250, res = 300, units = "px")
tmp <- summary_mat[,"kl_div"]; tmp <- tmp[tmp > 0]
xmin <- floor(log10(min(tmp)))
summary_mat[which(summary_mat[,"kl_div"] == 0),"kl_div"] <- 10^xmin
xmax <- ceiling(max(log10(summary_mat[,"kl_div"])))
plot(NA,
     xlim = c(xmin, xmax),
     ylim = c(0,1),
     bty = "n",
     xaxt = "n",
     xlab = "Separability (Median -Log10(p-value))",
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

points(log10(summary_mat[,"kl_div"]), summary_mat[,"r_squared"], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 1))
graphics.off()


########### ## ggrepel version of the plot
df <- as.data.frame(summary_mat)
df$Name <- rownames(df)
p1 <- ggplot2::ggplot(df, ggplot2::aes(x = kl_div, y = r_squared))
p1 <- p1 + ggplot2::geom_point()
p1 <- p1 + ggplot2::scale_x_log10()
# p1 <- p1 + ggplot2::xlim(1, max(df$kl_div))
p1 <- p1 + ggplot2::ylim(0, 1)
p1 <- p1 + ggrepel::geom_text_repel(data = df, ggplot2::aes(label = Name))
p1 <- p1 + Seurat::NoLegend()
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_adt_exploration2.png",
                p1, device = "png", width = 5, height = 5, units = "in")

###################################

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
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_dcca_adtUmap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

for(i in 1:length(antibody_list)){
  nrow_val <- ceiling(sqrt(length(antibody_list[[i]])))
  ncol_val <- ceiling(length(antibody_list[[i]])/nrow_val)
  plot1 <- Seurat::FeaturePlot(pbmc, 
                               reduction = "adt.umap",
                               features = antibody_list[[i]],
                               ncol = ncol_val)
  
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_dcca_adtUmap_", names(antibody_list)[i], "_antibodies.png"),
                  plot1, device = "png", width = 4*ncol_val, height = 4*nrow_val, units = "in")
}

for(i in 1:length(antibody_list)){
  df <- as.data.frame(summary_mat)
  df$Labeling <- rep(0, nrow(df))
  df$Labeling[which(rownames(df) %in% antibody_list[[i]])] <- 1
  df$Labeling <- as.factor(df$Labeling)
  df$Name <- rownames(df)
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = kl_div, y = r_squared))
  p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = Labeling))
  p1 <- p1 + ggplot2::scale_color_manual(values = c("black", "red"))
  p1 <- p1 + ggplot2::scale_x_log10()
  p1 <- p1 + ggplot2::ylim(0, 1)
  p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, Labeling == 1), 
                                      ggplot2::aes(label = Name, color = Labeling))
  p1 <- p1 + Seurat::NoLegend()
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_adt_exploration_", names(antibody_list)[i], ".png"),
                  p1, device = "png", width = 5, height = 5, units = "in")
}