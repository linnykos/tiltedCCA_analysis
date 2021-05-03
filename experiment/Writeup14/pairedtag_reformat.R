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


my_fcn2 <- function(dat2) {
  p <- progressr::progressor(along = 1:length(dat2))
  y <- purrr::reduce(dat2, function(x,y){
    p("x")
    cbind2(x,y)
  })
}
dat3 <- my_fcn2(dat2)
