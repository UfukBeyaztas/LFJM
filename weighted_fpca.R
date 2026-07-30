weighted_fpca <- function(
    functional_data,
    weights,
    proportion_variance = 0.95,
    number_components = NULL,
    center = TRUE
) {
  functional_data <- as.matrix(functional_data)
  weights <- as.numeric(weights)
  
  center_curve <- if (center) {
    colMeans(functional_data, na.rm = TRUE)
  } else {
    rep(0, ncol(functional_data))
  }
  
  centered_data <- sweep(
    functional_data,
    2,
    center_curve,
    "-"
  )
  centered_data[!is.finite(centered_data)] <- 0
  
  weighted_data <- sweep(
    centered_data,
    2,
    sqrt(weights),
    "*"
  )
  
  decomposition <- svd(weighted_data)
  eigenvalues <-
    decomposition$d^2 /
    max(1, nrow(functional_data) - 1)
  
  largest_eigenvalue <- if (length(eigenvalues) > 0) {
    max(eigenvalues)
  } else {
    0
  }
  
  if (!is.finite(largest_eigenvalue) ||
      largest_eigenvalue <= 1e-10) {
    return(
      list(
        scores = matrix(0, nrow(functional_data), 0),
        eigenfunctions = matrix(0, ncol(functional_data), 0),
        eigenvalues = numeric(0),
        center = center_curve
      )
    )
  }
  
  retained <- eigenvalues > largest_eigenvalue * 1e-10
  eigenvalues <- eigenvalues[retained]
  right_vectors <- decomposition$v[, retained, drop = FALSE]
  
  if (is.null(number_components)) {
    number_components <-
      which(
        cumsum(eigenvalues) / sum(eigenvalues) >=
          proportion_variance
      )[1]
  }
  
  number_components <- min(
    max(1, number_components),
    ncol(right_vectors)
  )
  
  right_vectors <- right_vectors[
    ,
    seq_len(number_components),
    drop = FALSE
  ]
  
  list(
    scores = weighted_data %*% right_vectors,
    eigenfunctions = sweep(
      right_vectors,
      1,
      sqrt(weights),
      "/"
    ),
    eigenvalues = eigenvalues[seq_len(number_components)],
    center = center_curve
  )
}