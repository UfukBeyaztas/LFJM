solve_system <- function(A, b, ridge = 1e-8) {
  A <- as.matrix(A)
  
  solution <- tryCatch(
    solve(A, b),
    error = function(e) NULL
  )
  
  if (!is.null(solution) && all(is.finite(solution))) {
    return(solution)
  }
  
  solution <- tryCatch(
    solve(A + ridge * diag(ncol(A)), b),
    error = function(e) NULL
  )
  
  if (!is.null(solution) && all(is.finite(solution))) {
    return(solution)
  }
  
  ginv(A) %*% b
}
