rm(list=ls())
genotype_est <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_est.rds")
genotype_values <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_values_h1_h2_c20.rds")
seg_table_filtered <- readRDS("../../../../data/Chiyun_SNU601/SNU601/seg_table_filtered.rds")
mat <- Matrix::readMM("../../../../data/Chiyun_SNU601/SNU601/matrix.mtx")
features <- read.csv("../../../../data/Chiyun_SNU601/SNU601/features.txt", header = F)
barcodes <- read.csv("../../../../data/Chiyun_SNU601/SNU601/barcodes.txt", header = F)
clone_assign <- readRDS("../../../../data/Chiyun_SNU601/SNU601/cloneAssign.rds")

rownames(mat) <- features[,1]
colnames(mat) <- barcodes[,1]
mat <- mat[,names(clone_assign)]


# see how many cells not assigned to any clones
length(which(is.na(clone_assign)))
table(clone_assign)

SNU <- Seurat::CreateSeuratObject(counts = mat,
                                  assay = "ATAC")
SNU[["clone"]] <- clone_assign
set.seed(10)
SNU <- Signac::RunTFIDF(SNU)
SNU <- Signac::FindTopFeatures(SNU, min.cutoff = 'q10')
SNU <- Seurat::ScaleData(SNU, vars.to.regress = "nCount_ATAC",
                         do.scale = F, do.center = F)
save(SNU, file = "../../../../out/Writeup14e/Writeup14e_SNU_atac_exploratory.RData")

SNU[["lsi"]] <- Signac::RunSVD(SNU[["ATAC"]]@scale.data)
set.seed(10)
SNU <- Seurat::RunUMAP(SNU, reduction = 'lsi', dims = 1:50, 
                       reduction.name = "umap.atac", 
                       reduction.key = "atacUMAP_")
save(SNU, file = "../../../../out/Writeup14e/Writeup14e_SNU_atac_exploratory.RData")

## now plot

# group_vec <- c("clone")
# for(group in group_vec){
#   plot1 <- Seurat::DimPlot(SNU, reduction = "umap.atac",
#                            group.by = group, label = TRUE,
#                            repel = TRUE, label.size = 2.5)
#   plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (ATAC)\n",group))
#   plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
#   ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_atac_exploration_", group, ".png"),
#                   plot1, device = "png", width = 5, height = 5, units = "in")
#   
# }
# 
# group_vec <- c("nCount_ATAC", "nFeature_ATAC")
# for(group in group_vec){
#   plot1 <- Seurat::FeaturePlot(SNU, reduction = "umap.atac",
#                            features = group)
#   plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (ATAC)\n",group))
#   plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
#   ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_atac_exploration_", group, ".png"),
#                   plot1, device = "png", width = 5, height = 5, units = "in")
#   
# }


