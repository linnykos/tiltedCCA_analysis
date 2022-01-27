rm(list=ls())
# load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca.RData")
# load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca_tmp.RData")
# dcca_res2 <- res
# class(dcca_res2) <- "dcca"
# 
# library(Seurat); library(Signac)
# library(multiomicCCA)
# 
# dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)
# 
# # now investigate celltypes
# print("Working on RNA graph")
# set.seed(10)
# rna_frnn <- multiomicCCA::construct_frnn(dcca_res2, 
#                                          data_1 = T, 
#                                          data_2 = F,
#                                          normalization_type = "cosine_itself")
# save(rna_frnn, 
#      file = "../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_clisi.RData")
# print(dim(rna_frnn$c_g))
# print(rna_frnn$c_g[1:5,1:5])
# 
# print("Working on ADT graph")
# set.seed(10)
# adt_frnn <- multiomicCCA::construct_frnn(dcca_res2, 
#                                          data_1 = F, 
#                                          data_2 = T,
#                                          normalization_type = "cosine_itself")
# save(rna_frnn, adt_frnn, 
#      file = "../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_clisi.RData")
# 
# print("Working on Common graph")
# set.seed(10)
# common_g <- multiomicCCA::combine_frnn(dcca_res2, 
#                                        g_1 = rna_frnn$c_g,
#                                        g_2 = adt_frnn$c_g,
#                                        nn = 30,
#                                        verbose = 1)
# save(rna_frnn, adt_frnn, common_g,
#      file = "../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_clisi.RData")
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_clisi.RData")

print("Working on RNA clisi")
set.seed(10)
rna_clisi <- multiomicCCA:::clisi_information(common_g, rna_frnn$d_g,
                                             membership_vec = factor(pbmc$celltype.l2),
                                             max_subsample_clisi = 4000)
rna_clisi$common_clisi$clisi_mat
rna_clisi$distinct_clisi$clisi_mat

save(rna_clisi, 
     file = "../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_clisi2.RData")

print("Working on ADT clisi")
set.seed(10)
adt_clisi <- multiomicCCA:::clisi_information(common_g, adt_frnn$d_g,
                                             membership_vec = factor(pbmc$celltype.l2),
                                             max_subsample_clisi = 4000)
adt_clisi$common_clisi$clisi_mat
adt_clisi$distinct_clisi$clisi_mat

save(rna_clisi, adt_clisi,
     file = "../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_clisi2.RData")
