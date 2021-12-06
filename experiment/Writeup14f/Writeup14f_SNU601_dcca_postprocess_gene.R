rm(list=ls())
load("../../../../out/Writeup14f/Writeup14f_SNU601_dcca.RData")
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2)

library(mclust); library(Seurat); library(Signac)

summary_mat <- compute_variable_summary(mat = dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2,
                                        common_mat = dcca_decomp$common_mat_2,
                                        metacell_clustering = factor(SNU$clone),
                                        verbose = 1)

# load in the NA mask
genotype_values <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_values_h1_h2_c20.rds")

mat_rho_org <- genotype_values[,grep("^rho*", colnames(genotype_values))]
mat_rho_org <- mat_rho_org[rownames(mat_rho_org) %in% rownames(SNU@meta.data),]
mat_theta_org <- genotype_values[,grep("^theta*", colnames(genotype_values))]
mat_theta_org <- mat_theta_org[rownames(mat_theta_org) %in% rownames(SNU@meta.data),]
all(rownames(mat_rho_org) == rownames(SNU@meta.data))
length(which(is.na(mat_rho_org)))
length(which(is.na(mat_theta_org)))
cna_org <- t(cbind(mat_rho_org, mat_theta_org))

# Wilcox rank-sum test
uniq_celltype <- sort(unique(SNU$clone))
combn_mat <- utils::combn(length(uniq_celltype), 2)
set.seed(10)
Seurat::Idents(SNU) <- "clone"
Seurat::DefaultAssay(SNU) <- "cna"
de_list <- lapply(1:ncol(combn_mat), function(i){
    ident_1 <- which(SNU$clone == uniq_celltype[combn_mat[1,i]])
    ident_2 <- which(SNU$clone == uniq_celltype[combn_mat[2,i]])
    
    p_val_vec <- sapply(1:nrow(SNU[["cna"]]@data), function(j){
        mask1 <- which(!is.na(cna_org[j,ident_1]))
        mask2 <- which(!is.na(cna_org[j,ident_2]))
        if(length(mask1) == 0 | length(mask2) == 0) return(1)
        
        vec1 <- SNU[["cna"]]@data[j,ident_1[mask1]]
        vec2 <- SNU[["cna"]]@data[j,ident_2[mask2]]
        
        tmp <- stats::wilcox.test(vec1, vec2)
        tmp$p.value
    })
    
    data.frame(variable = rownames(SNU[["cna"]]@data), 
               p_val = p_val_vec)
})
sapply(de_list, function(x){quantile(x[,2])})

summary_mat <- cbind(summary_mat, rep(0, nrow(summary_mat)))
colnames(summary_mat)[3] <- "p_val"
# overwrite column of summary_mat
for(i in 1:nrow(summary_mat)){
    gene_name <- rownames(summary_mat)[i]
    val <- mean(sapply(1:length(de_list), function(j){
        idx <- which(de_list[[j]][,1] == gene_name)
        if(length(idx) == 0) return(1)
        de_list[[j]][idx, "p_val"]
    }))
    
    summary_mat[i,"p_val"] <- val
}

# add jitter
set.seed(10)
summary_mat[,"p_val"] <- pmax(pmin(summary_mat[,"p_val"] + runif(nrow(summary_mat), min = -.05, max = .05), 1), 0)

col_vec <- rep(1, ncol(dcca_decomp$common_mat_2))
col_vec[grep("theta", colnames(dcca_decomp$common_mat_2))] <- 2

png("../../../../out/figures/Writeup14f/Writeup14f_SNU601_cna_exploration.png", 
    height = 1500, width = 1500, res = 300, units = "px")
plot(summary_mat[,3], summary_mat[,1], pch = 16,
     col = col_vec, xlab = "Separability (Average p-value, jittered)",
     ylab = "Alignment w/ common space (R^2)",
     main = "SNU601 (CNA):\n(Theta: red, Rho: black)")
graphics.off()

grep("-3:",  colnames(dcca_decomp$common_mat_2))
grep("-20:",  colnames(dcca_decomp$common_mat_2))
idx_3_rho <- c(7,8,9)
idx_3_theta <- c(71,72,73)
idx_20_rho <- c(59,60)
idx_20_theta <- c(123,124)

png("../../../../out/figures/Writeup14f/Writeup14f_SNU601_cna_exploration_specific.png", 
    height = 1500, width = 1500, res = 300, units = "px")
plot(summary_mat[,3], summary_mat[,1], pch = 16,
     col = "gray", xlab = "Separability (Average p-value, jittered)",
     ylab = "Alignment w/ common space (R^2)",
     main = "Chr3: red, Chr20: green\nRho: circle, Theta: triangle")
points(summary_mat[c(idx_3_rho, idx_3_theta, idx_20_rho, idx_20_theta),3], 
       summary_mat[c(idx_3_rho, idx_3_theta, idx_20_rho, idx_20_theta),1],
       pch = c(rep(16,3),rep(17,2), rep(16,3),rep(17,2)),
       col = c(rep(2, 5), rep(3, 5)),
       cex = 2)
graphics.off()

rownames(summary_mat)[order(summary_mat[,2], decreasing = T)[1:10]]

#####################
load("../../../../out/Writeup14f/Writeup14f_SNU_preprocessed.RData")
Seurat::DefaultAssay(SNU) <- "cna"
mat_2 <- Matrix::t(SNU[["cna"]]@data)

j_vec <- c(8,60)
uniq_clone <- sort(unique(SNU$clone))
color_vec <- scales::hue_pal()(length(uniq_clone))

for(j in 1:2){
    png(paste0("../../../../out/figures/Writeup14f/Writeup14f_SNU601_cna_exploration_specific_",
        rownames(summary_mat)[j_vec[j]], ".png"), 
        height = 2000, width = 1500, 
        res = 300, units = "px")
    par(mfrow = c(3,2), mar = c(4,4,4,0.5))
    xlim <- c(min(mat_2[,j_vec[j]]), quantile(mat_2[,j_vec[j]], prob = 0.95))
    for(i in 1:length(uniq_clone)){
        cell_idx <- intersect(
            intersect(which(SNU$clone == uniq_clone[i]),
                      which(mat_2[,j_vec[j]] >= xlim[1])),
            which(mat_2[,j_vec[j]] <= xlim[2]))
        hist(mat_2[cell_idx,j_vec[j]], col = color_vec[i], xlim = xlim,
             main = paste0("Clone ", uniq_clone[i]), breaks = 25)
    }
    graphics.off()
}
