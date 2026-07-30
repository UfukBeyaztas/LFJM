dual_basis_map <- function(eigenfunctions, grid) {
  weights <- trapezoidal_weights(grid)
  eigenfunctions <- as.matrix(eigenfunctions)
  
  if (ncol(eigenfunctions) == 0) {
    return(matrix(0, nrow(eigenfunctions), 0))
  }
  
  eigenfunctions %*%
    solve_system(
      crossprod(eigenfunctions, weights * eigenfunctions),
      diag(ncol(eigenfunctions))
    )
}