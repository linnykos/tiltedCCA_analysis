rm(list=ls())
library(Seurat); library(Signac)

load("../../../out/main/10x_greenleaf_lineageAlignment_pseudotime.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

##############

num_orderings <- 100
random_range <- 3
starting_celltype <- "Cyc. Prog."
verbose <- 1

rna_pseudotime <- greenleaf$rna_pseudotime
atac_pseudotime <- greenleaf$atac_pseudotime
celltype_vec <- greenleaf$celltype

cell_idx <- intersect(which(!is.na(rna_pseudotime)), which(!is.na(atac_pseudotime)))
cell_names_keep <- colnames(greenleaf)[cell_idx]

rna_pseudotime <- rna_pseudotime[cell_names_keep]
atac_pseudotime <- atac_pseudotime[cell_names_keep]
celltype_vec <- celltype_vec[cell_names_keep]

ordering_list <- lapply(1:num_orderings, function(i){
  if(verbose > 0) print(paste0("Working on ordering ", i))
  set.seed(10*i)
  
  starting_cell <- sample(which(celltype_vec == starting_celltype), 1)
  ordering_vec <- starting_cell
  names(ordering_vec) <- NULL
  
  while(TRUE){
    if(verbose > 1) print(paste0("Current length of ordering ", length(ordering_vec)))
    
    # find current pseudotime
    current_rna_pseudotime <- rna_pseudotime[ordering_vec[length(ordering_vec)]]
    current_atac_pseudotime <- atac_pseudotime[ordering_vec[length(ordering_vec)]]
    
    # find all later cells
    later_rna_idx <- which(rna_pseudotime > current_rna_pseudotime)
    later_atac_idx <- which(rna_pseudotime > current_rna_pseudotime)
    later_intersect_idx <- intersect(later_rna_idx, later_atac_idx)
    if(verbose > 1) print(paste0("Number of later cells ", length(later_intersect_idx)))
    if(length(later_intersect_idx) == 0) break()
    
    # construct the data frame
    df <- data.frame(idx = later_intersect_idx,
                     rna_pseudotime = rna_pseudotime[later_intersect_idx],
                     atac_pseudotime = atac_pseudotime[later_intersect_idx])
    df$rna_rank <- rank(df$rna_pseudotime)
    df$atac_rank <- rank(df$atac_pseudotime)
    df$max_rank <- pmax(df$rna_rank, df$atac_rank)
    
    min_idx <- df$idx[order(df$max_rank, decreasing = F)[1:random_range]]
    ordering_vec <- c(ordering_vec, sample(min_idx, 1))
    names(ordering_vec) <- NULL
  }
  
  if(verbose > 0) print(paste0("Ordering ", i, " has length ", length(ordering_vec)))
  ordering_vec
})

tmp <- rna_pseudotime[ordering_list[[30]]]; names(tmp) <- NULL; tmp

##########################

consensus_pseudotime <- greenleaf$consensus_pseudotime[cell_names_keep]
tcca_pseudotime <- greenleaf$tcca_pseudotime[cell_names_keep]

sequence_length_vec <- 2:10
num_trials <- 1000

.compute_alignment <- function(pseudotime_vec,
                               ordering_list, 
                               sequence_length_vec){
  vec <- sapply(1:length(sequence_length_vec), function(i){
    sequence_length <- sequence_length_vec[i]
    
    bool_vec <- sapply(1:num_trials, function(trial){
      set.seed(trial)
      
      ordering_idx <- sample(1:length(ordering_list), 1)
      selected_idx <- sort(sample(1:length(ordering_list[[ordering_idx]]), sequence_length))
      selected_cell_idx <- ordering_list[[ordering_idx]][selected_idx]
      
      rank_vec <- pseudotime_vec[selected_cell_idx]
      all(diff(rank_vec) >= 0)
    })
    
    length(which(bool_vec))/length(bool_vec)
  })
  
  names(vec) <- paste0("length_", sequence_length_vec)
  vec
}

consensus_score_vec <- .compute_alignment(pseudotime_vec = consensus_pseudotime,
                                          ordering_list = ordering_list, 
                                          sequence_length_vec = sequence_length_vec)

tcca_score_vec <- .compute_alignment(pseudotime_vec = tcca_pseudotime,
                                          ordering_list = ordering_list, 
                                          sequence_length_vec = sequence_length_vec)

cbind(consensus_score_vec, tcca_score_vec)