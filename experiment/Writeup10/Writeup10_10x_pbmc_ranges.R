rm(list=ls()); set.seed(10); gcinfo(TRUE)

library(Seurat); library(Signac); library(GenomeInfoDb)

## see LinkPeaks in https://github.com/timoast/signac/blob/master/R/links.R (line 198)

## https://github.com/timoast/signac/blob/master/R/links.R line 473
peak <- Signac::StringToGRanges(rownames(pbmc[["ATAC"]]), sep = c(":", "-"))
## https://github.com/timoast/signac/blob/master/R/utilities.R see GetTSSPositions line 438
gene_coords <- Signac:::CollapseToLongestTranscript(
  ranges = Signac::Annotation(object = pbmc[["ATAC"]])
)
genes <- Seurat::VariableFeatures(object = pbmc, assay = "SCT")
gene_coords_use <- gene_coords[gene_coords$gene_name %in% genes,]

dist_mat <- Signac:::DistanceToTSS(peak, gene_coords_use, distance = 500000)

# now let's quantify how many peaks are associatd with each gene
quantile(diff(dist_mat@p))
# we can also quantify how often which peaks are associated with genes
quantile(table(dist_mat@i))

# let's try this out for one gene, MS4A1
idx <- which(colnames(dist_mat) == "MS4A1")
peak_idx <- dist_mat@i[(dist_mat@p[idx]+1):(dist_mat@p[idx+1])]
peak[peak_idx]