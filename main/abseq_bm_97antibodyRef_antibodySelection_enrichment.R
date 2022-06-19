rm(list=ls())

library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/abseq_bm97Ref_tcca.RData")
load("../../../out/main/abseq_bm97Ref_varSelect_alternatives.RData")
bm_alt <- bm
load("../../../out/main/abseq_bm97Ref_varSelect.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

simplified_ct <- as.character(bm$ct)
simplified_ct[grep("CD4+", simplified_ct)] <- "CD4"
simplified_ct[grep("CD8+", simplified_ct)] <- "CD8"

celltype_vec <- c("CD4", "CD8")
# celltype_vec <- unique(simplified_ct)[c(grep("^CD4+",  unique(simplified_ct)), grep("^CD8+",  unique(simplified_ct)))]
# celltype_vec <- c("CD4+ memory T cells", "CD8+ central memory T cells", "CD8+ effector memory T cells", "CD8+CD103+ tissue resident memory T cells")
consensus_list <- vector("list", length = 4)
names(consensus_list) <- c("target", "alt_1", "alt_2", "alt_3")
consensus_list[[1]] <- consensus_pca$dimred_consensus[which(simplified_ct %in% celltype_vec),]
consensus_list[[2]] <- bm_alt[["consensusPCA1"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
consensus_list[[3]] <- bm_alt[["consensusPCA2"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
consensus_list[[4]] <- bm_alt[["consensusPCA3"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
membership_vec <- droplevels(as.factor(simplified_ct[which(simplified_ct %in% celltype_vec)]))

enrichment_list <- lapply(1:4, function(i){
  print(paste0("Working on ", i))
  consensus_embedding <- consensus_list[[i]]
  set.seed(10)
  enrichment_selected <- tiltedCCA:::postprocess_cell_enrichment(input_obj = consensus_embedding,
                                                                 membership_vec = membership_vec, 
                                                                 num_neigh = 100,
                                                                 bool_cosine = T,
                                                                 bool_intersect = F,
                                                                 max_subsample = 5000,
                                                                 min_deg = 50,
                                                                 verbose = 3)
})

for(i in 1:4){
  print(i)
  res <- enrichment_list[[i]]
  print(res$enrichment$df)
}

######################

celltype_vec <- c("Class switched memory B cells", "Nonswitched memory B cells", "Mature naive B cells")
consensus_list <- vector("list", length = 4)
names(consensus_list) <- c("target", "alt_1", "alt_2", "alt_3")
consensus_list[[1]] <- consensus_pca$dimred_consensus[which(simplified_ct %in% celltype_vec),]
consensus_list[[2]] <- bm_alt[["consensusPCA1"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
consensus_list[[3]] <- bm_alt[["consensusPCA2"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
consensus_list[[4]] <- bm_alt[["consensusPCA3"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
membership_vec <- droplevels(as.factor(simplified_ct[which(simplified_ct %in% celltype_vec)]))

enrichment_list <- lapply(1:4, function(i){
  print(paste0("Working on ", i))
  consensus_embedding <- consensus_list[[i]]
  set.seed(10)
  enrichment_selected <- tiltedCCA:::postprocess_cell_enrichment(input_obj = consensus_embedding,
                                                                 membership_vec = membership_vec, 
                                                                 num_neigh = 100,
                                                                 bool_cosine = T,
                                                                 bool_intersect = F,
                                                                 max_subsample = 5000,
                                                                 min_deg = 50,
                                                                 verbose = 3)
})

for(i in 1:4){
  print(i)
  res <- enrichment_list[[i]]
  print(res$enrichment$df)
}

