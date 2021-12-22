vec1 <- c(1, 0)
vec2 <- c(1/sqrt(2), 1/sqrt(2))

basis <- .representation_2d(vec1, vec2)
circle <- .construct_circle(basis$rep1, 
                            basis$rep2)
radian <- .compute_radian(circle = circle,
                          enforce_boundary = T,
                          percentage_val = 0.5, 
                          vec1 = basis$rep1,
                          vec2 = basis$rep2)
common_vec <- .position_from_circle(circle, radian)
# common_vec <- (1-tan(acos(vec1 %*% vec2)/2)) * (vec1 + vec2)/2
png("../../out/simulation/Writeup14e_simulation/example1_white.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(4,4,0.5,0.5))
plot_decomposition_2d(vec1, vec2, common_vec,
                      xlab = "2D plane, 1st dimension",
                      ylab = "2D plane, 2nd dimension",
                      plot_bg = F,
                      xlim = c(0,1), ylim = c(0,1), 
                      bty = "n",
                      length_arrow = 0.15,
                      lwd_arrow_black = 4,
                      lwd_arrow_redgreen = 4,
                      lwd_arrow_white = 5)
graphics.off()

png("../../out/simulation/Writeup14e_simulation/example1.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(4,4,0.5,0.5))
plot_decomposition_2d(vec1, vec2, common_vec,
                      xlab = "2D plane, 1st dimension",
                      ylab = "2D plane, 2nd dimension",
                      plot_bg = T,
                      bty = "n",
                      length_arrow = 0.15,
                      lwd_arrow_black = 4,
                      lwd_arrow_redgreen = 8,
                      lwd_arrow_white = 10)
graphics.off()


radian <- .compute_radian(circle = circle,
                          enforce_boundary = T,
                          percentage_val = 0.125, 
                          vec1 = basis$rep1,
                          vec2 = basis$rep2)
common_vec <- .position_from_circle(circle, radian)
# common_vec <- (1-tan(acos(vec1 %*% vec2)/2)) * (vec1 + vec2)/2
png("../../out/simulation/Writeup14e_simulation/example2.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(4,4,0.5,0.5))
plot_decomposition_2d(vec1, vec2, common_vec,
                      xlab = "2D plane, 1st dimension",
                      ylab = "2D plane, 2nd dimension",
                      plot_bg = T,
                      bty = "n",
                      length_arrow = 0.2,
                      lwd_arrow_black = 4,
                      lwd_arrow_redgreen = 8,
                      lwd_arrow_white = 10)
graphics.off()


radian <- .compute_radian(circle = circle,
                          enforce_boundary = T,
                          percentage_val = 1-0.125, 
                          vec1 = basis$rep1,
                          vec2 = basis$rep2)
common_vec <- .position_from_circle(circle, radian)
# common_vec <- (1-tan(acos(vec1 %*% vec2)/2)) * (vec1 + vec2)/2
png("../../out/simulation/Writeup14e_simulation/example3.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(4,4,0.5,0.5))
plot_decomposition_2d(vec1, vec2, common_vec,
                      xlab = "2D plane, 1st dimension",
                      ylab = "2D plane, 2nd dimension",
                      plot_bg = T,
                      bty = "n",
                      length_arrow = 0.15,
                      lwd_arrow_black = 4,
                      lwd_arrow_redgreen = 8,
                      lwd_arrow_white = 10)
graphics.off()


######################################

vec1 <- c(1, 0)
angle <- 75
vec2 <- c(cos(angle*pi/180), sin(angle*pi/180))
vec2 <- vec2/.l2norm(vec2)

common_vec <- (vec1+vec2)/2

png("../../out/simulation/Writeup14e_simulation/example_joint.png",
    height = 1200, width = 1200, units = "px", res = 300)
par(mar = c(4,4,0.5,0.5))
plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T,
     xlab = "2D plane, 1st dimension",
     ylab = "2D plane, 2nd dimension",
     bty = "n")

graphics::arrows(x0 = 0, y0 = 0, x1 = common_vec[1], y1 = common_vec[2], length = 0.1, col = "white", lwd = 3)
graphics::arrows(x1 = common_vec[1], y1 = common_vec[2], x0 = vec1[1], y0 = vec1[2], length = 0.1, col = "white", lwd = 3)
graphics::arrows(x1 = common_vec[1], y1 = common_vec[2], x0 = vec2[1], y0 = vec2[2], length = 0.1, col = "white", lwd = 3)

graphics::arrows(x0 = 0, y0 = 0, x1 = common_vec[1], y1 = common_vec[2], length = 0.1, col = 3, lwd = 3)
graphics::arrows(x1 = common_vec[1], y1 = common_vec[2], x0 = vec1[1], y0 = vec1[2], length = 0.1, col = 2, lwd = 3)
graphics::arrows(x1 = common_vec[1], y1 = common_vec[2], x0 = vec2[1], y0 = vec2[2], length = 0.1, col = 2, lwd = 3)

graphics::arrows(x0 = 0, y0 = 0, x1 = vec1[1], y1 = vec1[2], length = 0.1, lwd = 1.2)
graphics::arrows(x0 = 0, y0 = 0, x1 = vec2[1], y1 = vec2[2], length = 0.1, lwd = 1.2)
graphics.off()


vec1 <- c(1, 0)
vec2 <- c(1/sqrt(2), 1/sqrt(2))
basis <- .representation_2d(vec1, vec2)
circle <- .construct_circle(basis$rep1, 
                            basis$rep2)
radian <- .compute_radian(circle = circle,
                          enforce_boundary = T,
                          percentage_val = 0.5, 
                          vec1 = basis$rep1,
                          vec2 = basis$rep2)
common_vec <- .position_from_circle(circle, radian)

png("../../out/simulation/Writeup14e_simulation/example1_white_zoom.png",
    height = 700, width = 700, units = "px", res = 300)
par(mar = c(4,4,0.5,0.5))
plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T,
     xlab = "2D plane, 1st dim.",
     ylab = "2D plane, 2nd dim.",
     bty = "n")

graphics::arrows(x0 = 0, y0 = 0, x1 = common_vec[1], y1 = common_vec[2], length = 0.2, col = "white", lwd = 6)
graphics::arrows(x1 = common_vec[1], y1 = common_vec[2], x0 = vec1[1], y0 = vec1[2], length = 0.2, col = "white", lwd = 6)
graphics::arrows(x1 = common_vec[1], y1 = common_vec[2], x0 = vec2[1], y0 = vec2[2], length = 0.2, col = "white", lwd = 6)

graphics::arrows(x1 = common_vec[1], y1 = common_vec[2], x0 = vec1[1], y0 = vec1[2], length = 0.2, col = 2, lwd = 6)
graphics::arrows(x1 = common_vec[1], y1 = common_vec[2], x0 = vec2[1], y0 = vec2[2], length = 0.2, col = 2, lwd = 6)
graphics::arrows(x0 = 0, y0 = 0, x1 = common_vec[1], y1 = common_vec[2], length = 0.2, col = 3, lwd = 6)

graphics::arrows(x0 = 0, y0 = 0, x1 = vec1[1], y1 = vec1[2], length = 0.15, lwd = 3)
graphics::arrows(x0 = 0, y0 = 0, x1 = vec2[1], y1 = vec2[2], length = 0.15, lwd = 3)
graphics.off()


