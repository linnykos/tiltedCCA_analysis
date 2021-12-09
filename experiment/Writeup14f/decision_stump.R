decision_stump <- function(vec, label){
  stopifnot(is.factor(label), is.numeric(vec), 
            length(vec) == length(label),
            length(unique(label)) == 2)
  level_vec <- levels(label)
  
  # reorder
  n <- length(vec)
  idx <- order(vec)
  vec <- vec[idx]; label <- label[idx]
  
  # make into data frame
  pos_idx <- which(label == level_vec[1])
  pos_vec <- rep(0, n); pos_vec[pos_idx] <- 1
  neg_vec <- rep(1, n); neg_vec[-pos_idx] <- 0
  df <- data.frame(vec = vec, pos = pos_vec, neg = neg_vec)
  
  pos_cumsum <- cumsum(pos_vec)
  vec_cumsum <- cumsum(neg_vec)
}