predict.lfjm_lfpca <- function(
    object,
    functional_data,
    subject,
    time,
    ...
) {
  functional_data <- as.matrix(functional_data)
  grid <- object$grid
  number_grid_points <- ncol(functional_data)
  weights <- object$weights
  
  if (number_grid_points != length(grid)) {
    stop(
      "The new functional data must be evaluated on the fitted grid."
    )
  }
  
  if (nrow(functional_data) != length(subject) ||
      nrow(functional_data) != length(time)) {
    stop(
      "The rows of 'functional_data' must correspond to 'subject' and 'time'."
    )
  }
  
  if (object$presmooth) {
    functional_data <- t(
      apply(
        functional_data,
        1,
        function(curve) {
          smooth.spline(
            grid,
            curve,
            df = object$smooth_df
          )$y
        }
      )
    )
  }
  
  standardized_time <-
    (time - object$time_mean) /
    object$time_sd
  
  time_design <- polynomial_design(
    standardized_time,
    object$eta_degree
  )
  
  centered_data <-
    functional_data -
    time_design %*% object$eta_coef
  
  subject_ids <- sort(unique(subject))
  subject_index <- match(subject, subject_ids)
  number_subjects <- length(subject_ids)
  
  between_intercept <- matrix(
    0,
    number_subjects,
    number_grid_points
  )
  within_residuals <- matrix(
    0,
    nrow(functional_data),
    number_grid_points
  )
  
  if (object$use_slope) {
    between_slope <- matrix(
      0,
      number_subjects,
      number_grid_points
    )
    
    for (subject_number in seq_len(number_subjects)) {
      rows <- which(subject_index == subject_number)
      subject_curves <- centered_data[rows, , drop = FALSE]
      subject_time_design <-
        cbind(1, standardized_time[rows])
      
      subject_coefficients <- if (length(rows) >= 2) {
        ginv(
          crossprod(subject_time_design)
        ) %*%
          crossprod(
            subject_time_design,
            subject_curves
          )
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
        subject_time_design %*% subject_coefficients
    }
    
    stacked_between <- cbind(
      between_intercept,
      between_slope
    )
    stacked_center <- c(
      object$center_B0,
      object$center_B1
    )
    stacked_eigenfunctions <- rbind(
      object$phi0,
      object$phi1
    )
    stacked_weights <- c(weights, weights)
    
    between_scores <-
      (
        sweep(
          stacked_between,
          2,
          stacked_center,
          "-"
        ) *
          matrix(
            stacked_weights,
            number_subjects,
            2 * number_grid_points,
            byrow = TRUE
          )
      ) %*%
      stacked_eigenfunctions
  } else {
    for (subject_number in seq_len(number_subjects)) {
      rows <- which(subject_index == subject_number)
      subject_curves <- centered_data[rows, , drop = FALSE]
      subject_mean <- colMeans(subject_curves)
      
      between_intercept[subject_number, ] <- subject_mean
      within_residuals[rows, ] <-
        sweep(subject_curves, 2, subject_mean, "-")
    }
    
    between_scores <-
      (
        sweep(
          between_intercept,
          2,
          object$center_B0,
          "-"
        ) *
          matrix(
            weights,
            number_subjects,
            number_grid_points,
            byrow = TRUE
          )
      ) %*%
      object$phi0
  }
  
  within_scores <-
    (
      sweep(
        within_residuals,
        2,
        object$center_U,
        "-"
      ) *
        matrix(
          weights,
          nrow(within_residuals),
          number_grid_points,
          byrow = TRUE
        )
    ) %*%
    object$phiU
  
  rownames(between_scores) <- as.character(subject_ids)
  
  list(
    ids = subject_ids,
    subject_index = subject_index,
    scoreB_subject = between_scores,
    scoreU_visit = within_scores
  )
}