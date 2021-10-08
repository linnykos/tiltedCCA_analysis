rm(list=ls())
genotype_est <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_est.rds")
dim(genotype_est)
genotype_est[1:5,1:5]

genomat <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genomat.rds")
dim(genomat)
sort(unique(as.numeric(genomat)))
genomat[1:5,1:5]

seg_table_filtered <- readRDS("../../../../data/Chiyun_SNU601/SNU601/seg_table_filtered.rds")
dim(seg_table_filtered)
seg_table_filtered[1:5,]

mat <- Matrix::readMM("../../../../data/Chiyun_SNU601/SNU601/matrix.mtx")
dim(mat)
mat[1:5,1:5]
quantile(mat@x)
length(mat@x)/prod(dim(mat))

features <- read.csv("../../../../data/Chiyun_SNU601/SNU601/features.txt", header = F)
nrow(features)
features[1:5,1]

barcodes <- read.csv("../../../../data/Chiyun_SNU601/SNU601/barcodes.txt", header = F)
nrow(barcodes)
barcodes[1:5,1]

all(rownames(genomat) %in% barcodes[,1])
all(genomat[!is.na(genomat)] == genotype_est[!is.na(genotype_est)])
