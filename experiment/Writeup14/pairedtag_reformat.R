rm(list=ls())
library(tidyverse)
library(progressr)
progressr::handlers(global = TRUE)

dat <- readr::read_tsv("../../../../data/Pairedtag_mousebrain_RNA-Histone/exprMatrix.tsv")

# from https://stackoverflow.com/questions/52841870/how-to-cast-a-tibble-to-a-sparse-matrix
# dat2 <- purrr::map(dat[,-1], Matrix::Matrix, sparse = T) %>% purrr::reduce(cbind2)

# from https://cran.r-project.org/web/packages/progressr/vignettes/progressr-intro.html
my_fcn <- function(dat2) {
  p <- progressr::progressor(along = 1:ncol(dat2))
  y <- purrr::map(dat2, function(x){
    p("x")
    Matrix::Matrix(x, sparse = T)
  })
}

dat2 <- my_fcn(dat[,-1])

dat3 <- dat2[[1]]
vec_x <- unlist(pbapply::pblapply(dat2, function(x){x@x}))
vec_i <- unlist(pbapply::pblapply(dat2, function(x){x@i}))
vec_p <- c(0, cumsum(pbapply::pbsapply(dat2, function(x){x@p[2]})))
names(vec_x) <- NULL
names(vec_i) <- NULL
names(vec_p) <- NULL
vec_p <- as.integer(vec_p)
dat3@x <- vec_x; dat3@i <- vec_i; dat3@p <- vec_p; dat3@Dim[2] <- as.integer(length(dat2))
dat3@Dimnames[[1]] <- sapply(dat[,1], as.character)
dat3@Dimnames[[2]] <- colnames(dat)[-1]

# dat3[100:105,500:505]
# dat[100:105,501:506]

dat_tmp <- dat
dat <- dat3

save(dat, file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/dat.RData")
