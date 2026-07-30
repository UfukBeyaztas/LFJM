lfjm_control <- function(
    pve_between = 0.95,
    pve_within = 0.95,
    components_between = NULL,
    components_within = NULL,
    mean_degree = 3,
    presmooth = TRUE,
    smoothing_df = 12,
    use_functional_slope = FALSE,
    quadrature_points = 32,
    max_iterations = 500,
    multistart = TRUE,
    standard_errors = TRUE
) {
  if (!is.numeric(pve_between) || length(pve_between) != 1 ||
      !is.finite(pve_between) || pve_between <= 0 || pve_between > 1) {
    stop("'pve_between' must be in (0, 1].")
  }
  
  if (!is.numeric(pve_within) || length(pve_within) != 1 ||
      !is.finite(pve_within) || pve_within <= 0 || pve_within > 1) {
    stop("'pve_within' must be in (0, 1].")
  }
  
  if (!is.numeric(mean_degree) || length(mean_degree) != 1 ||
      !is.finite(mean_degree) || mean_degree < 0) {
    stop("'mean_degree' must be a non-negative number.")
  }
  
  if (!is.numeric(quadrature_points) || length(quadrature_points) != 1 ||
      !is.finite(quadrature_points) || quadrature_points < 2) {
    stop("'quadrature_points' must be at least 2.")
  }
  
  if (!is.numeric(max_iterations) || length(max_iterations) != 1 ||
      !is.finite(max_iterations) || max_iterations < 1) {
    stop("'max_iterations' must be positive.")
  }
  
  structure(
    list(
      pve_between = pve_between,
      pve_within = pve_within,
      components_between = components_between,
      components_within = components_within,
      mean_degree = mean_degree,
      presmooth = isTRUE(presmooth),
      smoothing_df = smoothing_df,
      use_functional_slope = isTRUE(use_functional_slope),
      quadrature_points = as.integer(quadrature_points),
      max_iterations = as.integer(max_iterations),
      multistart = isTRUE(multistart),
      standard_errors = isTRUE(standard_errors)
    ),
    class = "lfjm_control"
  )
}
