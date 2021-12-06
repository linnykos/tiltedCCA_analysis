rm(list=ls())
load("../../../../out/Writeup14f/Writeup14f_10x_mouseembryo_dcca.RData")
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)

library(mclust); library(Seurat); library(Signac)

summary_mat <- compute_variable_summary(mat = dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1,
                                        common_mat = dcca_decomp$common_mat_1,
                                        metacell_clustering = factor(mbrain2$label_Savercat),
                                        verbose = 2)

save.image(file = "../../../../out/Writeup14f/Writeup14f_10x_mouseembryo_dcca_gene.RData")

##############

load("../../../../out/Writeup14f/Writeup14f_10x_mouseembryo_dcca_gene.RData")
library(MAST); library(Seurat); library(Signac)
uniq_celltype <- sort(unique(mbrain2$label_Savercat))
label_vec <- c("Cortical", "Forebrain", "Glioblast", "Neuroblast", "Oligodendrocyte", "Radial")
for(i in 1:length(uniq_celltype)){
  newgroup <- rep(NA, nrow(mbrain2@meta.data))
  names(newgroup) <- rownames(mbrain2@meta.data)
  idx <- which(mbrain2$label_Savercat == uniq_celltype[i])
  newgroup[idx] <- paste0("in_", label_vec[i])
  newgroup[-idx] <- paste0("out_", label_vec[i])
  mbrain2[[paste0("indicator_", label_vec[i])]] <- newgroup
  Seurat::Idents(mbrain2) <- paste0("indicator_", label_vec[i])
}

combn_mat <- utils::combn(length(uniq_celltype), 2)
set.seed(10)
Seurat::Idents(mbrain2) <- "label_Savercat"
Seurat::DefaultAssay(mbrain2) <- "RNA"
de_list <- lapply(1:ncol(combn_mat), function(j){
  print(j)
  ident_1 <- uniq_celltype[combn_mat[1,j]]
  ident_2 <- uniq_celltype[combn_mat[2,j]]
  set.seed(10)
  Seurat::FindMarkers(mbrain2,
                      ident.1 = ident_1,
                      ident.2 = ident_2,
                      test.use = "MAST",
                      verbose = T)
})

save(uniq_celltype, combn_mat, de_list,
     file = "../../../../out/Writeup14f/Writeup14f_10x_mouseembryo_dcca_gene_de.RData")

hk_genes <- read.csv("../../../../data/housekeeping.txt")[,1]
rownames(summary_mat) <- toupper(rownames(summary_mat))
for(j in 1:length(de_list)){rownames(de_list[[j]]) <- toupper(rownames(de_list[[j]]))}
length(intersect(hk_genes, rownames(summary_mat)))
length(intersect(cc.genes$s.genes, rownames(summary_mat)))
length(intersect(cc.genes$g2m.genes, rownames(summary_mat)))

summary_mat <- cbind(summary_mat, rep(0, nrow(summary_mat)))
colnames(summary_mat)[3] <- "p_val"
# overwrite column of summary_mat
for(i in 1:nrow(summary_mat)){
  gene_name <- rownames(summary_mat)[i]
  val <- mean(sapply(1:length(de_list), function(j){
    idx <- which(rownames(de_list[[j]]) == gene_name)
    if(length(idx) == 0) return(1)
    de_list[[j]][idx, "p_val"]
  }))
  
  summary_mat[i,"p_val"] <- val
}

# add jitter
set.seed(10)
summary_mat[,"p_val"] <- pmax(pmin(summary_mat[,"p_val"] + runif(nrow(summary_mat), min = -.05, max = .05), 1), 0)

summary_mat[which(rownames(summary_mat) %in% hk_genes),]
summary_mat[which(rownames(summary_mat) %in% cc.genes$s.genes),]
summary_mat[which(rownames(summary_mat) %in% cc.genes$g2m.genes),]

negative_list <- list(Housekeeping = hk_genes, 
                      S = cc.genes$s.genes, 
                      G2M = cc.genes$g2m.genes)

# find the relevant indices
p_val_thres <- 1e-10
log_thres <- 5
num_instances <- 4
gene_list <- lapply(1:length(uniq_celltype), function(i){
  idx <- which(combn_mat == i, arr.ind = T)[,2]
  lis <- lapply(de_list[idx], function(mat){
    idx1 <- which(abs(mat[,2]) >= log_thres)
    idx2 <- which(abs(mat[,5]) <= p_val_thres)
    rownames(mat)[intersect(idx1, idx2)]
  })
  vec <- unlist(lis)
  tab <- table(vec)
  names(tab)[which(tab >= num_instances)]
})
sapply(gene_list, length)


png("../../../../out/figures/Writeup14f/Writeup14f_10x_mouseembryo_rna_exploration.png",
    height = 1500, width = 1500, res = 300, units = "px")
plot(summary_mat[,3], summary_mat[,1], pch = 16,
     col = rgb(0.5, 0.5, 0.5, 0.5), xlab = "Separability (Average p-value, jittered)",
     ylab = "Alignment w/ common space (R^2)",
     main = "Mouse embryo (RNA)")
graphics.off()


color_vec <- scales::hue_pal()(length(unique(mbrain2$label_Savercat)))
for(i in 1:length(uniq_celltype)){
  idx <- which(rownames(summary_mat) %in% gene_list[[i]])

  png(paste0("../../../../out/figures/Writeup14f/Writeup14f_10x_mouseembryo_rna_exploration_",
             label_vec[i], ".png"),
      height = 1500, width = 1500, res = 300, units = "px")
  plot(summary_mat[,3], summary_mat[,1], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 0.5), xlab = "Separability (Average p-value, jittered)",
       ylab = "Alignment w/ common space (R^2)",
       main = paste0("Mouse embryo (RNA): ", label_vec[i]))

  # quantile_vec <- quantile(summary_mat[idx,1], probs = c(0.25, 0.5, 0.75))
  quantile_vec <- mean(summary_mat[idx,1]) + c(-1,0,1)*sd(summary_mat[idx,1])

  polygon(cbind(c(-1e3, -1e3, 1e3, 1e3),
                c(quantile_vec[1], quantile_vec[3], quantile_vec[3], quantile_vec[1])),
          density = 10, col = color_vec[i])

  lines(c(-1e5, 1e5), rep(quantile_vec[2], 2), lty = 2, lwd = 2, col = color_vec[i])

  points(summary_mat[idx,3], summary_mat[idx,1], pch = 16,
         col = "white", cex = 2)
  points(summary_mat[idx,3], summary_mat[idx,1], pch = 16,
         col = color_vec[i], cex = 1.5)
  graphics.off()
}

for(i in 1:length(negative_list)){
  idx <- which(rownames(summary_mat) %in% negative_list[[i]])
  
  png(paste0("../../../../out/figures/Writeup14f/Writeup14f_10x_mouseembryo_rna_exploration_",
             names(negative_list)[i], ".png"),
      height = 1500, width = 1500, res = 300, units = "px")
  plot(summary_mat[,3], summary_mat[,1], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 0.5), xlab = "Separability (Average p-value, jittered)",
       ylab = "Alignment w/ common space (R^2)",
       main = paste0("Mouse embryo (RNA): ", names(negative_list)[i], " genes"))
  
  # quantile_vec <- quantile(summary_mat[idx,1], probs = c(0.25, 0.5, 0.75))
  quantile_vec <- mean(summary_mat[idx,1]) + c(-1,0,1)*sd(summary_mat[idx,1])
  
  polygon(cbind(c(-1e3, -1e3, 1e3, 1e3),
                c(quantile_vec[1], quantile_vec[3], quantile_vec[3], quantile_vec[1])),
          density = 10, col = 2)
  
  lines(c(-1e5, 1e5), rep(quantile_vec[2], 2), lty = 2, lwd = 2, col = 2)
  
  points(summary_mat[idx,3], summary_mat[idx,1], pch = 16,
         col = "white", cex = 2)
  points(summary_mat[idx,3], summary_mat[idx,1], pch = 16,
         col = 2, cex = 1.5)
  graphics.off()
}


