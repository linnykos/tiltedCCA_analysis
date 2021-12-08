rm(list=ls())
load("../../../../out/Writeup14f/Writeup14f_citeseq_bm_dcca.RData")
source("gene_exploration.R")
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)

library(Seurat); library(Signac)

tab <- table(bm$celltype.l2)
cell_lists <- list("CD14 Mono", "CD16 Mono", "CD4 Memory", "CD4 Naive",
                   "CD8 Effector_1", "CD8 Effector_2", "CD8 Memory_1", "CD8 Memory_2",
                   "CD8 Naive", c("MAIT", "gdT"), c("Memory B", "Naive B"),
                   "NK", "pDC", "Prog_RBC")
celltype.l3 <- rep(NA, nrow(dcca_decomp$common_mat_2))
names(celltype.l3) <- rownames(dcca_decomp$common_mat_2)
for(i in 1:length(cell_lists)){
  idx <- which(bm$celltype.l2 %in% cell_lists[[i]])
  celltype.l3[idx] <- i
}
bm$celltype.l3 <- celltype.l3
Seurat::Idents(bm) <- "celltype.l3"

summary_mat2 <- compute_variable_summary(mat = dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2,
                                        common_mat = dcca_decomp$common_mat_2,
                                        metacell_clustering = factor(bm$celltype.l2),
                                        verbose = 1)

combn_mat <- utils::combn(length(cell_lists), 2)
set.seed(10)
Seurat::DefaultAssay(bm) <- "ADT"
adt_de_list <- lapply(1:ncol(combn_mat), function(i){
  print(i)
  ident_1 <- which(bm$celltype.l3 == as.character(combn_mat[1,i]))
  ident_2 <- which(bm$celltype.l3 == as.character(combn_mat[2,i]))
  
  p_val_vec <- sapply(1:nrow(bm[["ADT"]]@data), function(j){
    vec1 <- bm[["ADT"]]@data[j,ident_1]
    vec2 <- bm[["ADT"]]@data[j,ident_2]
    
    tmp <- stats::wilcox.test(vec1, vec2)
    tmp$p.value
  })
  
  data.frame(variable = rownames(bm[["ADT"]]@data), 
             p_val = p_val_vec)
})
summary_mat2 <- cbind(summary_mat2, rep(0, nrow(summary_mat2)))
colnames(summary_mat2)[3] <- "p_val"
# overwrite column of summary_mat
for(i in 1:nrow(summary_mat2)){
  gene_name <- rownames(summary_mat2)[i]
  val <- mean(sapply(1:length(adt_de_list), function(j){
    idx <- which(adt_de_list[[j]][,1] == gene_name)
    if(length(idx) == 0) return(1)
    adt_de_list[[j]][idx, "p_val"]
  }))
  
  summary_mat2[i,"p_val"] <- val
}

summary_mat1 <- compute_variable_summary(mat = dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1,
                                         common_mat = dcca_decomp$common_mat_1,
                                         metacell_clustering = factor(bm$celltype.l2),
                                         verbose = 1)
set.seed(10)
Seurat::DefaultAssay(bm) <- "RNA"
Seurat::Idents(bm) <- "celltype.l3"
gene_de_list <- lapply(1:ncol(combn_mat), function(i){
  print(j)
  ident_1 <- as.character(combn_mat[1,i])
  ident_2 <- as.character(combn_mat[2,i])
  
  set.seed(10)
  Seurat::FindMarkers(bm,
                      ident.1 = ident_1,
                      ident.2 = ident_2,
                      test.use = "MAST",
                      verbose = F)
})

save(summary_mat1, summary_mat2, adt_de_list, gene_de_list,
     file = "../../../../out/Writeup14f/Writeup14f_citeseq_bm_dcca_postprocess_gene.RData")

summary_mat1 <- cbind(summary_mat1, rep(0, nrow(summary_mat1)))
colnames(summary_mat1)[3] <- "p_val"
# overwrite column of summary_mat
for(i in 1:nrow(summary_mat1)){
  gene_name <- rownames(summary_mat1)[i]
  val <- mean(sapply(1:length(gene_de_list), function(j){
    idx <- which(gene_de_list[[j]][,1] == gene_name)
    if(length(idx) == 0) return(1)
    gene_de_list[[j]][idx, "p_val"]
  }))
  
  summary_mat1[i,"p_val"] <- val
}

save(summary_mat1, summary_mat2, adt_de_list, gene_de_list,
     file = "../../../../out/Writeup14f/Writeup14f_citeseq_bm_dcca_postprocess_gene.RData")
