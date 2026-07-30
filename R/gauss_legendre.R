gauss_legendre <- function(number_nodes) {
  indices <- seq_len(number_nodes - 1)
  off_diagonal <- indices / sqrt(4 * indices^2 - 1)
  Jacobi_matrix <- matrix(0, number_nodes, number_nodes)
  
  Jacobi_matrix[cbind(indices, indices + 1)] <- off_diagonal
  Jacobi_matrix[cbind(indices + 1, indices)] <- off_diagonal
  
  eigen_decomposition <- eigen(Jacobi_matrix, symmetric = TRUE)
  ordering <- order(eigen_decomposition$values)
  
  list(
    nodes = (eigen_decomposition$values[ordering] + 1) / 2,
    weights = eigen_decomposition$vectors[1, ordering]^2
  )
}
