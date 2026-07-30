estimate_mean_surface <- function(
    functional_data,
    time,
    time_degree = 4
) {
  functional_data <- as.matrix(functional_data)
  time <- as.numeric(time)
  
  polynomial_degree <- min(
    3,
    max(0, length(unique(time)) - 1),
    time_degree
  )
  
  time_design <- polynomial_design(time, polynomial_degree)
  coefficients <- qr.coef(qr(time_design), functional_data)
  coefficients[!is.finite(coefficients)] <- 0
  
  list(
    fitted_values = time_design %*% coefficients,
    coefficients = coefficients,
    degree = polynomial_degree
  )
}
