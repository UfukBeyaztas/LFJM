lfjm_lfpca <- function(
    functional_data,
    subject,
    time,
    grid,
    pve_between = 0.95,
    pve_within = 0.95,
    components_between = NULL,
    components_within = NULL,
    mean_degree = 4,
    presmooth = TRUE,
    smoothing_df = 12,
    use_functional_slope = FALSE
) {
  functional_data <- as.matrix(functional_data)
  subject <- as.vector(subject)
  time <- as.numeric(time)
  grid <- as.numeric(grid)
  
  if (nrow(functional_data) != length(subject) ||
      nrow(functional_data) != length(time)) {
    stop(
      "The rows of 'functional_data' must correspond to 'subject' and 'time'."
    )
  }
  
  if (ncol(functional_data) != length(grid)) {
    stop(
      "The columns of 'functional_data' must correspond to the functional grid."
    )
  }
  
  subject_ids <- sort(unique(subject))
  subject_index <- match(subject, subject_ids)
  number_grid_points <- ncol(functional_data)
  
  if (presmooth) {
    functional_data <- t(
      apply(
        functional_data,
        1,
        function(curve) {
          smooth.spline(
            grid,
            curve,
            df = smoothing_df
          )$y
        }
      )
    )
  }
  
  time_mean <- mean(time)
  time_sd <- sd(time)
  
  if (!is.finite(time_sd) || time_sd < 1e-8) {
    time_sd <- 1
  }
  
  standardized_time <- (time - time_mean) / time_sd
  
  mean_surface <- estimate_mean_surface(
    functional_data,
    standardized_time,
    mean_degree
  )
  
  centered_data <-
    functional_data -
    mean_surface$fitted_values
  
  number_subjects <- length(subject_ids)
  between_intercept <- matrix(
    0,
    number_subjects,
    number_grid_points
  )
  between_slope <- matrix(
    0,
    number_subjects,
    number_grid_points
  )
  within_residuals <- matrix(
    0,
    nrow(functional_data),
    number_grid_points
  )
  
  for (subject_number in seq_len(number_subjects)) {
    rows <- which(subject_index == subject_number)
    subject_curves <- centered_data[rows, , drop = FALSE]
    
    if (use_functional_slope) {
      time_design <- cbind(1, standardized_time[rows])
      
      subject_coefficients <- if (length(rows) >= 2) {
        ginv(crossprod(time_design)) %*%
          crossprod(time_design, subject_curves)
      } else {
        rbind(
          subject_curves[1, ],
          rep(0, number_grid_points)
        )
      }
      
      between_intercept[subject_number, ] <-
        subject_coefficients[1, ]
      between_slope[subject_number, ] <-
        subject_coefficients[2, ]
      within_residuals[rows, ] <-
        subject_curves -
        time_design %*% subject_coefficients
    } else {
      subject_mean <- colMeans(subject_curves)
      between_intercept[subject_number, ] <- subject_mean
      within_residuals[rows, ] <-
        sweep(subject_curves, 2, subject_mean, "-")
    }
  }
  
  weights <- trapezoidal_weights(grid)
  
  between_fpca <- if (use_functional_slope) {
    weighted_fpca(
      cbind(between_intercept, between_slope),
      c(weights, weights),
      pve_between,
      components_between
    )
  } else {
    weighted_fpca(
      between_intercept,
      weights,
      pve_between,
      components_between
    )
  }
  
  within_fpca <- weighted_fpca(
    within_residuals,
    weights,
    pve_within,
    components_within
  )
  
  number_between_components <- ncol(between_fpca$scores)
  number_within_components <- ncol(within_fpca$scores)
  
  if (use_functional_slope) {
    phi0 <- between_fpca$eigenfunctions[
      seq_len(number_grid_points),
      ,
      drop = FALSE
    ]
    phi1 <- between_fpca$eigenfunctions[
      number_grid_points + seq_len(number_grid_points),
      ,
      drop = FALSE
    ]
    center_between_slope <- between_fpca$center[
      number_grid_points + seq_len(number_grid_points)
    ]
  } else {
    phi0 <- between_fpca$eigenfunctions
    phi1 <- matrix(
      0,
      number_grid_points,
      number_between_components
    )
    center_between_slope <- rep(0, number_grid_points)
  }
  
  phiU <- within_fpca$eigenfunctions
  center_between_intercept <-
    between_fpca$center[seq_len(number_grid_points)]
  
  reconstructed_between_intercept <-
    between_fpca$scores %*% t(phi0) +
    matrix(
      center_between_intercept,
      number_subjects,
      number_grid_points,
      byrow = TRUE
    )
  
  reconstructed_between_slope <-
    between_fpca$scores %*% t(phi1) +
    matrix(
      center_between_slope,
      number_subjects,
      number_grid_points,
      byrow = TRUE
    )
  
  reconstructed_within <-
    within_fpca$scores %*% t(phiU) +
    matrix(
      within_fpca$center,
      nrow(functional_data),
      number_grid_points,
      byrow = TRUE
    )
  
  reconstructed_between_visit <-
    reconstructed_between_intercept[
      subject_index,
      ,
      drop = FALSE
    ] +
    standardized_time *
    reconstructed_between_slope[
      subject_index,
      ,
      drop = FALSE
    ]
  
  rownames(between_fpca$scores) <-
    as.character(subject_ids)
  rownames(reconstructed_between_intercept) <-
    as.character(subject_ids)
  
  result <- list(
    call = match.call(),
    ids = subject_ids,
    subject_index = subject_index,
    grid = grid,
    sgrid = grid,
    weights = weights,
    phi0 = phi0,
    phi1 = phi1,
    phiU = phiU,
    scoreB_subject = between_fpca$scores,
    scoreU_visit = within_fpca$scores,
    B0_hat_subject = reconstructed_between_intercept,
    B_hat_visit = reconstructed_between_visit,
    U_hat_visit = reconstructed_within,
    eigB = between_fpca$eigenvalues,
    eigU = within_fpca$eigenvalues,
    KB = number_between_components,
    KU = number_within_components,
    time_mean = time_mean,
    time_sd = time_sd,
    eta_coef = mean_surface$coefficients,
    eta_degree = mean_surface$degree,
    center_B0 = center_between_intercept,
    center_B1 = center_between_slope,
    center_U = within_fpca$center,
    presmooth = presmooth,
    smooth_df = smoothing_df,
    use_slope = use_functional_slope
  )
  
  class(result) <- "lfjm_lfpca"
  result
}