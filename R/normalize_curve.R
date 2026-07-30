normalize_curve <- function(curve, weights) {
  curve_norm <- sqrt(sum(weights * curve^2))
  
  if (!is.finite(curve_norm) || curve_norm < 1e-12) {
    return(curve)
  }
  
  curve / curve_norm
}
