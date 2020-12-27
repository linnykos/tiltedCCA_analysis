# try 1) separate PCAs, 2) concatenate both matrices,
# 3) concatenate the principal loadings, 4) the vanilla WNN,
# 5) the vanilla D-CCA

rm(list=ls())
set.seed(10)

library(Seurat); library(Signac)
load("../../out/Writeup9_10x_pbmc_dimred_all.RData")

#############################
# first define cell-types
# from https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html
pbmc <- Seurat::FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- Seurat::RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- Seurat::FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE) ## see http://www.ludowaltman.nl/slm/ for SLM (smart local moving)
pbmc <- Seurat::FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)
Seurat::Idents(pbmc) <- "sub.cluster"

pbmc <- Seurat::RenameIdents(pbmc, '19' = 'pDC','20' = 'HSPC','15' = 'cDC')
pbmc <- Seurat::RenameIdents(pbmc, '0' = 'CD14 Mono', '8' ='CD14 Mono', '5' = 'CD16 Mono')
pbmc <- Seurat::RenameIdents(pbmc, '17' = 'Naive B', '11' = 'Intermediate B', '10' = 'Memory B', '21' = 'Plasma')
pbmc <- Seurat::RenameIdents(pbmc, '7' = 'NK')
pbmc <- Seurat::RenameIdents(pbmc, '4' = 'CD4 TCM', '13'= "CD4 TEM", '3' = "CD4 TCM", '16' ="Treg", '1' ="CD4 Naive", '14' = "CD4 Naive")
pbmc <- Seurat::RenameIdents(pbmc, '2' = 'CD8 Naive', '9'= "CD8 Naive", '12' = 'CD8 TEM_1', '6_0' = 'CD8 TEM_2', '6_1' ='CD8 TEM_2')
pbmc <- Seurat::RenameIdents(pbmc, '18' = 'MAIT')
pbmc <- Seurat::RenameIdents(pbmc, '6_2' ='gdT', '6_3' = 'gdT')
pbmc$celltype <- Seurat::Idents(pbmc)

head(pbmc@meta.data)

##########################################3

# first try separate PCAs, already calculated
png("../../out/Writeup10_pbmc_10x_seurat_rna_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "celltype", 
                 label = TRUE, label.size = 2.5, repel = TRUE) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA UMAP (Vanilla)")
graphics.off()

png("../../out/Writeup10_pbmc_10x_seurat_atac_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "celltype", 
                 label = TRUE, label.size = 2.5, repel = TRUE) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("ATAC UMAP (Vanilla)")
graphics.off()

png("../../out/Writeup10_pbmc_10x_seurat_wnn_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", 
                 label = TRUE, label.size = 2.5, repel = TRUE) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("WNN UMAP (Seurat)")
graphics.off()

##################################

# check if scale.data really is scaled (it is)
mat <- pbmc[["ATAC"]]@scale.data
mat[70:80,90:100]
mean(mat[1,]); sd(mat[1,]); length(which(abs(mat[1,]) == 10))
zz <- which(apply(mat, 1, function(x){any(x != 0) & length(which(abs(x)==10)) == 0}))
head(zz)
mean(mat[8,]); sd(mat[8,]); length(which(abs(mat[8,]) == 10))
class(mat)
rm(list=c("mat"))

# pick a rank
mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ATAC"]]@scale.data)
dim(mat_1); dim(mat_2)
svd_res_1 <- RSpectra::svds(mat_1, k = 25)
svd_res_2 <- RSpectra::svds(mat_2, k = 25)
-diff(svd_res_1$d)/svd_res_1$d[-1]; -diff(svd_res_2$d)/svd_res_2$d[-1]
K <- 7
png("../../out/Writeup10_pbmc_10x_rank.png", height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,2))
graphics::plot(svd_res_1$d, xlab = "Rank", ylab = "Singular value", main = "SVD for RNA", pch = 16)
graphics::lines(rep(K+.5, 2), c(-1e5, 1e5), col = "red", lty = 2, lwd = 2)
graphics::plot(svd_res_2$d, xlab = "Rank", ylab = "Singular value", main = "SVD for ATAC", pch = 16)
graphics::lines(rep(K+.5, 2), c(-1e5, 1e5), col = "red", lty = 2, lwd = 2)
graphics.off()

##############################
# try method 2, just concenating the matrices naively
mat_total <- cbind(mat_1/svd_res_1$d[1], mat_2/svd_res_2$d[1]) # oops! error, "cannot allocate vector of size 8.3 Gb"

