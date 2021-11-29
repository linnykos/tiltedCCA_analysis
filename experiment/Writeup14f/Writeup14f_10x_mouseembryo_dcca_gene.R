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
library(MAST)
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

set.seed(10)
Seurat::DefaultAssay(mbrain2) <- "SCT"
de_list <- lapply(1:length(uniq_celltype), function(i){
  Seurat::Idents(mbrain2) <- paste0("indicator_", label_vec[i])
  ident_1 <- paste0("in_", label_vec[i])
  ident_2 <- paste0("out_", label_vec[i])
  set.seed(10)
  Seurat::FindMarkers(mbrain2,
                      slot = "scale.data",
                      ident.1 = ident_1,
                      ident.2 = ident_2,
                      test.use = "wilcox",
                      verbose = T)
})

png("../../../../out/figures/Writeup14f/Writeup14f_10x_mouseembryo_rna_exploration.png", 
    height = 1500, width = 1500, res = 300, units = "px")
plot(summary_mat[,2], summary_mat[,1], pch = 16,
     col = rgb(0.5, 0.5, 0.5, 0.5), xlab = "Separability (KL divergence)",
     ylab = "Alignment w/ common space (R^2)",
     main = "Mouse embryo (RNA)")
graphics.off()


color_vec <- scales::hue_pal()(length(unique(mbrain2$label_Savercat)))
for(i in 1:length(uniq_celltype)){
  gene_names <- rownames(de_list[[i]])[1:50]
  idx <- which(rownames(summary_mat) %in% gene_names)
  
  png(paste0("../../../../out/figures/Writeup14f/Writeup14f_10x_mouseembryo_rna_exploration_",
             label_vec[i], ".png"), 
      height = 1500, width = 1500, res = 300, units = "px")
  plot(summary_mat[,2], summary_mat[,1], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 0.5), xlab = "Separability (KL divergence)",
       ylab = "Alignment w/ common space (R^2)",
       main = "Mouse embryo (RNA)")
  points(summary_mat[idx,2], summary_mat[idx,1], pch = 16,
         col = "white", cex = 2)
  points(summary_mat[idx,2], summary_mat[idx,1], pch = 16,
         col = color_vec[i], cex = 1.5)
  graphics.off()
}


