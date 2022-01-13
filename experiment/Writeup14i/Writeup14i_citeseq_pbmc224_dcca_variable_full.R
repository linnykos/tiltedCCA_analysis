rm(list=ls())
load("../../../../out/Writeup14i/Writeup14i_citeseq_pbmc224_dcca.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

library(MAST); library(Seurat); library(Signac)
uniq_celltype <- as.character(sort(unique(pbmc$celltype.l2)))
combn_mat <- utils::combn(length(uniq_celltype), 2)

set.seed(10)
Seurat::Idents(pbmc) <- "celltype.l2"
Seurat::DefaultAssay(pbmc) <- "ADT"
de_list <- lapply(1:ncol(combn_mat), function(j){
  print(paste0(j, " of ", ncol(combn_mat)))
  
  ident_1 <- uniq_celltype[combn_mat[1,j]]
  ident_2 <- uniq_celltype[combn_mat[2,j]]
  set.seed(10)
  Seurat::FindMarkers(pbmc,
                      features = rownames(pbmc[["ADT"]]),
                      ident.1 = ident_1,
                      ident.2 = ident_2,
                      test.use = "wilcox",
                      min.pct = 0,
                      logfc.threshold = 0,
                      only.pos = F,
                      verbose = T)
})

save(uniq_celltype, combn_mat, de_list,
     file = "../../../../out/Writeup14i/Writeup14i_citeseq_pbmc224_dcca_variable_full_de.RData")
