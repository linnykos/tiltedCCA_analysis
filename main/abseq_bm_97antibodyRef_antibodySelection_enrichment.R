rm(list=ls())

library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/abseq_bm97Ref_tcca.RData")
load("../../../out/main/abseq_bm97Ref_varSelect_alternatives.RData")
bm_alt <- bm
load("../../../out/main/abseq_bm97Ref_varSelect.RData")
source("bm_97antibodyRef_colorPalette.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

simplified_ct <- as.character(bm$ct)

# celltype_vec <- unique(simplified_ct)[grep("^CD8+", unique(simplified_ct))]
celltype_vec <- c("CD8+ naive T cells", "CD8+ central memory T cells", "CD8+CD103+ tissue resident memory T cells")

consensus_list <- vector("list", length = 5)
names(consensus_list) <- c("target", "alt_1", "alt_2", "rna", "original")
consensus_list[[1]] <- consensus_pca$dimred_consensus[which(simplified_ct %in% celltype_vec),]
consensus_list[[2]] <- bm_alt[["consensusPCA1"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
consensus_list[[3]] <- bm_alt[["consensusPCA2"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
consensus_list[[4]] <- bm[["pca"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),1:20]
consensus_list[[5]] <- bm[["consensusPCA"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
membership_vec <- droplevels(as.factor(simplified_ct[which(simplified_ct %in% celltype_vec)]))
k <- length(consensus_list)

enrichment_list <- lapply(1:k, function(i){
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

for(i in 1:k){
  print(i)
  res <- enrichment_list[[i]]
  print(res$enrichment$df)
}

#################################

celltype_vec <- unique(simplified_ct)[grep("^CD4+", unique(simplified_ct))]

consensus_list <- vector("list", length = 5)
names(consensus_list) <- c("target", "alt_1", "alt_2", "rna", "original")
consensus_list[[1]] <- consensus_pca$dimred_consensus[which(simplified_ct %in% celltype_vec),]
consensus_list[[2]] <- bm_alt[["consensusPCA1"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
consensus_list[[3]] <- bm_alt[["consensusPCA2"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
consensus_list[[4]] <- bm[["pca"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),1:20]
consensus_list[[5]] <- bm[["consensusPCA"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
membership_vec <- droplevels(as.factor(simplified_ct[which(simplified_ct %in% celltype_vec)]))
k <- length(consensus_list)

enrichment_list <- lapply(1:k, function(i){
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

for(i in 1:k){
  print(i)
  res <- enrichment_list[[i]]
  print(res$enrichment$df)
}

main <- "CD4"
for(i in 1:3){
  vec <- enrichment_list[[i]]$enrichment$df[,"value"]/enrichment_list[[5]]$enrichment$df[,"value"]
  names(vec) <- enrichment_list[[i]]$enrichment$df[,"celltype"]
  
  png(paste0("../../../out/figures/main/abseq_bm97Ref_varSelect_histogram_", main, "_", names(consensus_list)[[i]], ".png"),
      height = 930, width = 560, units = "px", res = 500)
  par(mar = c(0.5,1,0.5,0))
  barplot(vec,  ylim = c(0, 1), space = 0,
          col = col_palette[names(vec)], 
          names.arg = rep("", length(vec)),
          xaxt = "n", yaxt = "n", bty = "n")
  graphics::axis(2,
                 cex.axis = 2,
                 lwd = 2,
                 lwd.ticks = 2)
  graphics.off()
}

######################

celltype_vec <- c("Class switched memory B cells", "Nonswitched memory B cells", "Mature naive B cells", "Immature B cells")
consensus_list <- vector("list", length = 5)
names(consensus_list) <- c("target", "alt_1", "alt_2", "rna", "original")
consensus_list[[1]] <- consensus_pca$dimred_consensus[which(simplified_ct %in% celltype_vec),]
consensus_list[[2]] <- bm_alt[["consensusPCA1"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
consensus_list[[3]] <- bm_alt[["consensusPCA2"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
consensus_list[[4]] <- bm[["pca"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),1:20]
consensus_list[[5]] <- bm[["consensusPCA"]]@cell.embeddings[which(simplified_ct %in% celltype_vec),]
membership_vec <- droplevels(as.factor(simplified_ct[which(simplified_ct %in% celltype_vec)]))
k <- length(consensus_list)

enrichment_list <- lapply(1:k, function(i){
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

for(i in 1:k){
  print(i)
  res <- enrichment_list[[i]]
  print(res$enrichment$df)
}


main <- "B"
for(i in 1:3){
  vec <- enrichment_list[[i]]$enrichment$df[,"value"]/enrichment_list[[5]]$enrichment$df[,"value"]
  names(vec) <- enrichment_list[[i]]$enrichment$df[,"celltype"]
  
  png(paste0("../../../out/figures/main/abseq_bm97Ref_varSelect_histogram_", main, "_", names(consensus_list)[[i]], ".png"),
      height = 930, width = 560, units = "px", res = 500)
  par(mar = c(0.5,1,0.5,0))
  barplot(vec,  ylim = c(0, 1), space = 0,
          col = col_palette[names(vec)], 
          names.arg = rep("", length(vec)),
          xaxt = "n", yaxt = "n", bty = "n")
  graphics::axis(2,
                 cex.axis = 2,
                 lwd = 2,
                 lwd.ticks = 2)
  graphics.off()
}


