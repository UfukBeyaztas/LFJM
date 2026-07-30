orthogonalize_columns <- function(A, weights) {
  A <- as.matrix(A)
  
  for (k in seq_len(ncol(A))) {
    column <- A[, k]
    
    if (k > 1) {
      for (l in seq_len(k - 1)) {
        column <- column -
          sum(weights * column * A[, l]) * A[, l]
      }
    }
    
    A[, k] <- normalize_curve(column, weights)
  }
  
  A
}