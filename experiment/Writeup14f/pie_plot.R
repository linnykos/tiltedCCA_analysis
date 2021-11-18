pie_custom <- function(x, offset = c(0,0), edges = 200, radius = 0.8, 
                       clockwise = T, 
                       init.angle = 90, 
                       col = 1:length(x), 
                       border = c(rep(2, length(x)/2), rep(1, length(x)/2)), lwd = 1){
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  if (length(border) == 1) border <- rep_len(border, nx)
  if (!is.null(lwd)) lwd <- rep_len(lwd, nx)
  twopi <- ifelse(clockwise, -2 * pi, 2 * pi)
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    graphics::polygon(c(P$x, 0) + offset[1], c(P$y, 0) + offset[2], 
                      border = border[i], col = col[i], lwd = lwd[i])
    P <- t2xy(mean(x[i + 0:1]))
  }
  
  invisible()
}