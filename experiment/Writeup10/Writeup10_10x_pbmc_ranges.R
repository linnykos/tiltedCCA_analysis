rm(list=ls()); set.seed(10); gcinfo(TRUE)

library(Seurat); library(Signac); library(GenomeInfoDb)
load("../../out/Writeup10_10x_pbmc_preprocess4.RData")
load("../../out/Writeup10_10x_pbmc_dcca_metacells2.RData")

## see LinkPeaks in https://github.com/timoast/signac/blob/master/R/links.R (line 198)

## https://github.com/timoast/signac/blob/master/R/links.R line 473
peak <- Signac::StringToGRanges(rownames(pbmc[["ATAC"]]), sep = c(":", "-"))
## https://github.com/timoast/signac/blob/master/R/utilities.R see GetTSSPositions line 438
gene_coords <- Signac:::CollapseToLongestTranscript(
  ranges = Signac::Annotation(object = pbmc[["ATAC"]])
)
genes <- Seurat::VariableFeatures(object = pbmc, assay = "SCT")
gene_coords_use <- gene_coords[gene_coords$gene_name %in% genes,]

link_mat <- Signac:::DistanceToTSS(peak, gene_coords_use, distance = 500000)

dim(link_mat)

# now let's quantify how many peaks are associatd with each gene
quantile(diff(link_mat@p))
# we can also quantify how often which peaks are associated with genes
quantile(table(link_mat@i))

# let's try this out for one gene, MS4A1
idx <- which(colnames(link_mat) == "MS4A1")
peak_idx <- link_mat@i[(link_mat@p[idx]+1):(link_mat@p[idx+1])]
peak[peak_idx]

length(link_mat@i)

##################

# for each gene-peak combination, compute two things:
# 1) the correlation b/w gene and peak in the common space, 
# 2) the combined explained variability 

# first, get everything in check
res <- dcca_decomposition(dcca_res, rank_c = K)
colnames(res$common_mat_1) <- Seurat::VariableFeatures(object = pbmc, assay = "SCT")
colnames(res$distinct_mat_1) <- Seurat::VariableFeatures(object = pbmc, assay = "SCT")
colnames(res$common_mat_2) <- rownames(pbmc[["ATAC"]]@scale.data)
colnames(res$distinct_mat_2) <- rownames(pbmc[["ATAC"]]@scale.data)

vec1 <- colnames(res$common_mat_1); vec2 <- colnames(link_mat)
subset_idx <- which(vec1 %in% vec2)
vec1 <- vec1[subset_idx]
gene_idx <- subset_idx[.matching_idx(vec1, vec2)]
gene_tab <- data.frame(gene_name = colnames(link_mat), idx = gene_idx)
head(gene_tab); head(colnames(res$common_mat_1)[gene_tab$idx])

head(rownames(link_mat)); head(colnames(res$common_mat_2))
vec1 <- colnames(res$common_mat_2); vec2 <- rownames(link_mat)
peak_idx <- .matching_idx(vec1, vec2)
peak_tab <- data.frame(peak_name = rownames(link_mat), idx = peak_idx)
head(peak_tab); head(colnames(res$common_mat_2)[peak_tab$idx])

len <- length(link_mat@i)

summary_mat <- sapply(1:len, function(x){
  print(x)
  
  residual <- x - link_mat@p
  idx <- which(residual > 0)
  i <- which.min(residual[idx]) 
  j <- link_mat@i[x]
  val1 <- stats::cor(res$common_mat_1[,gene_tab[i,2]], res$common_mat_2[,peak_tab[i,2]], method = "spearman")
  
  c1 <- .l2norm(res$common_mat_1[,gene_tab[i,2]]); d1 <- .l2norm(res$distinct_mat_1[,gene_tab[i,2]])
  c2 <- .l2norm(res$common_mat_2[,peak_tab[i,2]]); d2 <- .l2norm(res$distinct_mat_2[,peak_tab[i,2]])
  val2 <- mean(c(c1/(c1+d1), c2/(c2+d2)))
  
  c(val1, val2)
}) # about 15 minutes

summary_mat <- t(summary_mat)

quantile(summary_mat[,1]); quantile(summary_mat[,2])

vec1 <- summary_mat[,1]
vec2 <- summary_mat[,2]
vec2 <- sign(vec1)*vec2
vec1 <- abs(vec1)^2

png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_linkage.png", height = 1500, width = 1500, units = "px", res = 300)
plot(vec2, vec1, xlab = "% variability explained", ylab = "Spearman correlation squared",
     main = "Gene-peak linkage in common space", 
     ylim = c(0,1),
     pch = 16, col = rgb(0.5, 0.5, 0.5, 0.2))
lines(rep(0,2), c(-10,10), col = "red", lty = 2)
graphics.off()


# save.image("../../out/Writeup10_10x_pbmc_dcca_summary.RData")

