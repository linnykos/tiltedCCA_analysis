rm(list=ls())
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")
load("../../../out/main/citeseq_pbmc224_varSelect_alternatives.RData")
pbmc_alt <- pbmc
load("../../../out/main/citeseq_pbmc224_varSelect.RData")


library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

membership_vec <- droplevels(as.factor(pbmc$celltype.l2))
print("Enrichment selected")
set.seed(10)
enrichment_selected <- postprocess_cell_enrichment(input_obj = consensus_pca$dimred_consensus,
                                                   membership_vec = membership_vec, 
                                                   num_neigh = multiSVD_obj$param$snn_num_neigh,
                                                   bool_cosine = multiSVD_obj$param$snn_bool_cosine,
                                                   bool_intersect = multiSVD_obj$param$snn_bool_intersect,
                                                   max_subsample = 1000,
                                                   min_deg = multiSVD_obj$param$snn_min_deg,
                                                   verbose = 2)

print("Enrichment alternative 1")
set.seed(10)
enrichment_alt_1 <- postprocess_cell_enrichment(input_obj = pbmc_alt[["adt.umap1"]],
                                                membership_vec = membership_vec, 
                                                num_neigh = multiSVD_obj$param$snn_num_neigh,
                                                bool_cosine = multiSVD_obj$param$snn_bool_cosine,
                                                bool_intersect = multiSVD_obj$param$snn_bool_intersect,
                                                max_subsample = 1000,
                                                min_deg = multiSVD_obj$param$snn_min_deg,
                                                verbose = 2)

print("Enrichment alternative 2")
set.seed(10)
enrichment_alt_2 <- postprocess_cell_enrichment(input_obj = pbmc_alt[["adt.umap2"]],
                                                membership_vec = membership_vec, 
                                                num_neigh = multiSVD_obj$param$snn_num_neigh,
                                                bool_cosine = multiSVD_obj$param$snn_bool_cosine,
                                                bool_intersect = multiSVD_obj$param$snn_bool_intersect,
                                                max_subsample = 1000,
                                                min_deg = multiSVD_obj$param$snn_min_deg,
                                                verbose = 2)

print("Enrichment alternative 3")
set.seed(10)
enrichment_alt_3 <- postprocess_cell_enrichment(input_obj = pbmc_alt[["adt.umap3"]],
                                                membership_vec = membership_vec, 
                                                num_neigh = multiSVD_obj$param$snn_num_neigh,
                                                bool_cosine = multiSVD_obj$param$snn_bool_cosine,
                                                bool_intersect = multiSVD_obj$param$snn_bool_intersect,
                                                max_subsample = 1000,
                                                min_deg = multiSVD_obj$param$snn_min_deg,
                                                verbose = 2)

save(enrichment_selected, enrichment_alt_1,
     enrichment_alt_2, enrichment_alt_3,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224_varSelect_enrichment.RData")


