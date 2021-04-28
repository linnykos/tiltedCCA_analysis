rm(list=ls())

set.seed(10)
n_clust <- 100
B_mat1 <- matrix(c(0.9, 0, 0, 
                   0, 0.9, 0,
                   0, 0, 0.9), 3, 3, byrow = T)
B_mat2 <- matrix(c(0.9, 0.85, 0, 
                   0.85, 0.9, 0,
                   0, 0, 1), 3, 3, byrow = T)
K <- ncol(B_mat1)
membership_vec <- c(rep(1, n_clust), rep(2, n_clust), rep(3, n_clust))
n <- length(membership_vec); true_membership_vec <- membership_vec
svd_u_1 <- multiomicCCA::generate_sbm_orthogonal(B_mat1, membership_vec, centered = T)
svd_u_2 <- multiomicCCA::generate_sbm_orthogonal(B_mat2, membership_vec, centered = T)

p_1 <- 20; p_2 <- 40
svd_d_1 <- sqrt(n*p_1)*c(1.5,1); svd_d_2 <- sqrt(n*p_2)*c(1.5,1)
svd_v_1 <- multiomicCCA::generate_random_orthogonal(p_1, K-1)
svd_v_2 <- multiomicCCA::generate_random_orthogonal(p_2, K-1)

dat <- multiomicCCA::generate_data(svd_u_1, svd_u_2, svd_d_1, svd_d_2, svd_v_1, svd_v_2, 
                                   noise_val = 0.1)
K <- 2
dcca_res <- dcca_factor(dat$mat_1, dat$mat_2, rank_1 = K, rank_2 = K, verbose = F)
dcca_decomp <- dcca_decomposition(dcca_res, rank_c = K, verbose = F)

clisi_1 <- clisi_information(dcca_decomp$common_mat_1, dcca_decomp$distinct_mat_1,
                         membership_vec = as.factor(membership_vec),
                         rank_c = K, rank_d = K, nn = 50, frnn_approx = 0, 
                         min_subsample_cell = n_clust,
                         verbose = F)
clisi_2 <- clisi_information(dcca_decomp$common_mat_2, dcca_decomp$distinct_mat_2,
                                membership_vec = as.factor(membership_vec),
                                rank_c = K, rank_d = K, nn = 50, frnn_approx = 0,
                             min_subsample_cell = n_clust, 
                                verbose = F)

########

col_vec = scales::hue_pal()(nrow(clisi_1$common_clisi$membership_info))
par_mar = c(4,2.5,0.5,0.5)
par_oma = c(0,0,2,0)
asp = T
pch_main = 16
cex_main = 1.5
pch_bg = 16
cex_bg = 1
alpha_bg = 0.5
l_bg = 95
c_bg = 50
xlim = c(0,1)
ylim = c(0,1)
gridsize = 5
col_grid = grDevices::rgb(0.8,0.8,0.8)
lty_grid = 3
lwd_grid = 1
col_diag = "firebrick"
lty_diag = 2
lwd_diag = 2
xlab1 = "Distinct information 1"
xlab2 = "Distinct information 2"
ylab = "Common information"
ylab_dist = 0.5
main = "cLISI Information"
cex_text_main = 1.5

stopifnot(class(clisi_1) == "clisi", class(clisi_2) == "clisi",
          all(dim(clisi_1$common_clisi$cell_info) == dim(clisi_2$common_clisi$cell_info)),
          all(dim(clisi_1$common_clisi$membership_info) == dim(clisi_2$common_clisi$membership_info)),
          all(clisi_1$common_clisi$cell_info$celltype == clisi_2$common_clisi$cell_info$celltype))

bg_col_vec <- scales::alpha(scales::col2hcl(col_vec, l = l_bg, c = c_bg), alpha = alpha_bg)
graphics::par(mfrow = c(1,2), mar = par_mar, oma = par_oma)
x_vec <- seq(xlim[1], xlim[2], length.out = gridsize)
y_vec <- seq(ylim[1], ylim[2], length.out = gridsize)

graphics::plot(NA, xlim = sort(-1*xlim), ylim = ylim, asp = asp, xlab = xlab1,
               ylab = "", yaxt = 'n', xaxt = 'n')
graphics::title(ylab = ylab, line = ylab_dist)
graphics::axis(1, at = -1*x_vec, labels = as.character(round(x_vec, 2)))
graphics::axis(4, at = y_vec, labels = NA)
.draw_grid(x_vec, y_vec, xlim, ylim, col_grid, lty_grid, lwd_grid,
           flip = T)

.plot_clisi_cell(clisi_1, bg_col_vec, pch_bg, cex_bg, flip = T)
graphics::lines(c(0,-1),c(0,1), col = col_diag, lty = lty_diag, lwd = lwd_diag)
.plot_clisi_type(clisi_1, col_vec, pch_main, cex_main, flip = T)

####

graphics::plot(NA, xlim = xlim, ylim = ylim, asp = asp, xlab = xlab2,
               ylab = "", yaxt = 'n', xaxt = 'n')
graphics::axis(1, at = x_vec, labels = as.character(round(x_vec, 2)))
graphics::axis(2, at = y_vec, labels = as.character(round(y_vec, 2)))
.draw_grid(x_vec, y_vec, xlim, ylim, col_grid, lty_grid, lwd_grid,
           flip = F)

.plot_clisi_cell(clisi_2, bg_col_vec, pch_bg, cex_bg, flip = F)
graphics::lines(c(0,1),c(0,1), col = col_diag, lty = lty_diag, lwd = lwd_diag)
.plot_clisi_type(clisi_2, col_vec, pch_main, cex_main, flip = F)

####

graphics::mtext(main, outer = TRUE, cex = cex_text_main)