rm(list=ls())
load("../../../../out/Writeup14f/Writeup14f_10x_mouseembryo_dcca.RData")
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)

library(mclust); library(Seurat); library(Signac)

summary_mat <- compute_variable_summary(mat = dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1,
                                        common_mat = dcca_decomp$common_mat_1,
                                        metacell_clustering = factor(mbrain2$label_Savercat),
                                        verbose = 2)

save.image(file = "../../../../out/Writeup14f/Writeup14f_10x_mouseembryo_dcca_gene.RData")