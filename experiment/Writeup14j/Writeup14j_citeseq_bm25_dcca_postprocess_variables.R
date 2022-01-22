rm(list=ls())
load("../../../../out/Writeup14f/Writeup14f_citeseq_bm_dcca_postprocess_gene.RData")
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

# fix the p-values in summary_mat2
for(i in 1:nrow(summary_mat2)){
  gene_name <- rownames(summary_mat2)[i]
  val <- quantile(sapply(1:length(adt_de_list), function(j){
    idx <- which(adt_de_list[[j]][,1] == gene_name)
    if(length(idx) == 0) return(1)
    adt_de_list[[j]][idx, "p_val"]
  }), probs = 0.5)
  
  summary_mat2[i,"p_val"] <- val
}

# fix the p-values in summary_mat1
for(i in 1:nrow(summary_mat1)){
  if(i %% floor(nrow(summary_mat1)/10) == 0) cat('*')
  
  gene_name <- rownames(summary_mat1)[i]
  val <- quantile(sapply(1:length(gene_de_list), function(j){
    idx <- which(rownames(gene_de_list[[j]]) == gene_name)
    if(length(idx) == 0) return(1)
    gene_de_list[[j]][idx, "p_val"]
  }), probs = 0.5)
  
  summary_mat1[i,"p_val"] <- val
}

summary_mat1[,"p_val"] <- -log10(summary_mat1[,"p_val"])
summary_mat2[,"p_val"] <- -log10(summary_mat2[,"p_val"])

#########################

load("../../../../out/Writeup14j/Writeup14j_citeseq_bm25_dcca.RData")
source("../Writeup14f/gene_exploration.R")
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)

summary_mat1b <- compute_variable_summary(mat = mat_1b, 
                                          common_mat = dcca_decomp$common_mat_1,
                                          metacell_clustering = NA)
summary_mat2b <- compute_variable_summary(mat = mat_2b, 
                                          common_mat = dcca_decomp$common_mat_2,
                                          metacell_clustering = NA)
summary_mat1[,1] <- summary_mat1b[,1]
summary_mat2[,1] <- summary_mat2b[,1]

#################################
#################################
#################################

png("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm_adt_exploration.png",
    height = 1250, width = 1250, res = 300, units = "px")
tmp <- summary_mat2[,3]; tmp <- tmp[tmp > 0]
xmin <- floor(log10(min(tmp)))
summary_mat2[which(summary_mat2[,3] == 0),3] <- 10^xmin
xmax <- ceiling(max(log10(summary_mat2[,3])))
plot(NA,
     xlim = c(xmin, xmax),
     ylim = c(0,1),
     bty = "n",
     xaxt = "n",
     xlab = "Separability (Median -Log10(p-value))",
     ylab = "Alignment w/ common space (R^2)",
     main = "CITE-Seq: Bone marrow (Protein)")
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

points(log10(summary_mat2[,3]), summary_mat2[,1], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 1))
graphics.off()

########### ## ggrepel version of the plot
df <- as.data.frame(summary_mat2)
df$Name <- rownames(df)
p1 <- ggplot2::ggplot(df, ggplot2::aes(x = p_val, y = r_squared))
p1 <- p1 + ggplot2::geom_point()
p1 <- p1 + ggplot2::scale_x_log10()
p1 <- p1 + ggplot2::xlim(1, max(df$p_val))
p1 <- p1 + ggplot2::ylim(0, 1)
p1 <- p1 + ggrepel::geom_text_repel(data = df, ggplot2::aes(label = Name))
p1 <- p1 + Seurat::NoLegend()
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm_adt_exploration2.png",
                p1, device = "png", width = 5, height = 5, units = "in")

###########################

png("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm_rna_exploration.png",
    height = 1250, width = 1250, res = 300, units = "px")
tmp <- summary_mat1[,3]; tmp <- tmp[tmp > 0]
xmin <- floor(log10(min(tmp)))
summary_mat1[which(summary_mat1[,3] == 0),3] <- 10^xmin
xmax <- ceiling(max(log10(summary_mat1[,3])))
plot(NA,
     xlim = c(xmin, xmax),
     ylim = c(0,1),
     bty = "n",
     xaxt = "n",
     xlab = "Separability (Median -Log10(p-value))",
     ylab = "Alignment w/ common space (R^2)",
     main = "CITE-Seq: Bone marrow (RNA)")
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

points(log10(summary_mat1[,3]), summary_mat1[,1], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 0.2))
graphics.off()

############################################
cell_lists <- list("CD14 Mono", "CD16 Mono", "CD4 Memory", "CD4 Naive",
                   "CD8 Effector_1", "CD8 Effector_2", "CD8 Memory_1", "CD8 Memory_2",
                   "CD8 Naive", c("MAIT", "gdT"), c("Memory B", "Naive B"),
                   "NK", "pDC", "Prog_RBC")
combn_mat <- utils::combn(length(cell_lists), 2)

gene_de_list <- lapply(gene_de_list, function(mat){
  mat[rownames(mat) %in% rownames(summary_mat1),]
})
p_val_thres <- 1e-5
log_thres <- 3
num_instances <- round(length(cell_lists)/2)
gene_list <- lapply(1:length(cell_lists), function(i){
  idx <- which(combn_mat == i, arr.ind = T)[,2]
  lis <- lapply(gene_de_list[idx], function(mat){
    idx1 <- which(abs(mat[,2]) >= log_thres)
    idx2 <- which(abs(mat[,5]) <= p_val_thres)
    rownames(mat)[intersect(idx1, idx2)]
  })
  vec <- unlist(lis)
  tab <- table(vec)
  name_vec <- names(tab)[which(tab >= num_instances)]
  
  if(length(name_vec) <= 5){
    # find the genes that have the smallest median p-value
    all_genes <- unique(unlist(lapply(gene_de_list[idx], rownames)))
    all_genes <- all_genes[!all_genes %in% name_vec]
    
    p_val_vec <- sapply(all_genes, function(gene_name){
      val <- quantile(sapply(idx, function(j){
        row_idx <- which(rownames(gene_de_list[[j]]) == gene_name)
        if(length(row_idx) == 0) return(1)
        gene_de_list[[j]][row_idx, "p_val"]
      }), probs = 0.5)
    })
    
    name_vec <- c(name_vec, all_genes[order(p_val_vec, decreasing = F)][1:(5-length(name_vec))])
  }
  
  name_vec
})
sapply(gene_list, length)

hk_genes <- read.csv("../../../../data/housekeeping.txt")[,1]
negative_list <- list(Housekeeping = hk_genes[which(hk_genes %in% rownames(summary_mat1))], 
                      S = cc.genes$s.genes[which(cc.genes$s.genes %in% rownames(summary_mat1))], 
                      G2M = cc.genes$g2m.genes[which(cc.genes$g2m.genes %in% rownames(summary_mat1))])


label_vec <- sapply(cell_lists, function(x){paste0(x, collapse = "+")})
file_vec <- sapply(cell_lists, function(x){
  gsub(" ", "", paste0(x, collapse = "_"))
})
color_vec_tmp <- scales::hue_pal()(length(unique(bm$celltype.l2)))
color_vec <- sapply(cell_lists, function(x){
  color_vec_tmp[which(sort(unique(bm$celltype.l2)) == x[1])]
})
for(i in 1:length(cell_lists)){
  idx <- which(rownames(summary_mat1) %in% gene_list[[i]])
  
  png(paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm_rna_exploration_", 
             file_vec[i], ".png"),
      height = 1250, width = 1250, res = 300, units = "px")
  tmp <- summary_mat1[,3]; tmp <- tmp[tmp > 0]
  xmin <- floor(log10(min(tmp)))
  summary_mat1[which(summary_mat1[,3] == 0),3] <- 10^xmin
  xmax <- ceiling(max(log10(summary_mat1[,3])))
  plot(NA,
       xlim = c(xmin, xmax),
       ylim = c(0,1),
       bty = "n",
       xaxt = "n",
       xlab = "Separability (Median -Log10(p-value))",
       ylab = "Alignment w/ common space (R^2)",
       main = paste0("CITE-Seq: Bone marrow (RNA)\n", label_vec[i], " genes"))
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
  
  quantile_vec <- mean(summary_mat1[idx,1]) + c(-1,0,1)*sd(summary_mat1[idx,1])
  
  polygon(cbind(c(-1e3, -1e3, 1e3, 1e3),
                c(quantile_vec[1], quantile_vec[3], quantile_vec[3], quantile_vec[1])),
          density = 20, col = color_vec[i])
  
  lines(c(-1e5, 1e5), rep(quantile_vec[2], 2), lty = 2, lwd = 2, col = color_vec[i])
  
  points(log10(summary_mat1[,3]), summary_mat1[,1], pch = 16,
         col = rgb(0.5, 0.5, 0.5, 0.2))
  points(log10(summary_mat1[idx,3]), summary_mat1[idx,1], pch = 16,
         col = "white", cex = 2)
  points(log10(summary_mat1[idx,3]), summary_mat1[idx,1], pch = 16,
         col = color_vec[i], cex = 1.5)
  graphics.off()
}

for(i in 1:length(negative_list)){
  idx <- which(rownames(summary_mat1) %in% negative_list[[i]])
  
  png(paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_bm_rna_exploration_", 
             names(negative_list)[i], ".png"),
      height = 1250, width = 1250, res = 300, units = "px")
  tmp <- summary_mat1[,3]; tmp <- tmp[tmp > 0]
  xmin <- floor(log10(min(tmp)))
  summary_mat1[which(summary_mat1[,3] == 0),3] <- 10^xmin
  xmax <- ceiling(max(log10(summary_mat1[,3])))
  plot(NA,
       xlim = c(xmin, xmax),
       ylim = c(0,1),
       bty = "n",
       xaxt = "n",
       xlab = "Separability (Median -Log10(p-value))",
       ylab = "Alignment w/ common space (R^2)",
       main = paste0("CITE-Seq: Bone marrow (RNA)\n", names(negative_list)[i], " genes"))
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
  
  quantile_vec <- mean(summary_mat1[idx,1]) + c(-1,0,1)*sd(summary_mat1[idx,1])
  
  polygon(cbind(c(-1e3, -1e3, 1e3, 1e3),
                c(quantile_vec[1], quantile_vec[3], quantile_vec[3], quantile_vec[1])),
          density = 20, col = 2)
  
  lines(c(-1e5, 1e5), rep(quantile_vec[2], 2), lty = 2, lwd = 2, col = 2)
  
  points(log10(summary_mat1[,3]), summary_mat1[,1], pch = 16,
         col = rgb(0.5, 0.5, 0.5, 0.2))
  points(log10(summary_mat1[idx,3]), summary_mat1[idx,1], pch = 16,
         col = "white", cex = 2)
  points(log10(summary_mat1[idx,3]), summary_mat1[idx,1], pch = 16,
         col = 2, cex = 1.5)
  graphics.off()
}


