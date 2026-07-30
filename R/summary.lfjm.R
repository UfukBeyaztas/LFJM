summary.lfjm <- function(object, ...) {
  number_longitudinal_parameters <-
    length(object$beta)
  
  estimates <- c(
    object$beta,
    sigma_e = object$sigma_e,
    sigma_u = object$sigma_u,
    gamma0 = object$gamma0,
    gamma_Q = object$gamma_Q,
    alpha = object$alpha,
    Tau = object$Tau
  )
  
  standard_errors <- rep(
    NA_real_,
    length(estimates)
  )
  
  if (!is.null(object$covariance_matrix)) {
    raw_standard_errors <- sqrt(
      pmax(
        diag(object$covariance_matrix),
        0
      )
    )
    
    standard_errors[
      seq_len(number_longitudinal_parameters)
    ] <-
      raw_standard_errors[
        seq_len(number_longitudinal_parameters)
      ]
    
    standard_errors[
      number_longitudinal_parameters + 1
    ] <-
      object$sigma_e *
      raw_standard_errors[
        number_longitudinal_parameters + 1
      ]
    
    standard_errors[
      number_longitudinal_parameters + 2
    ] <-
      object$sigma_u *
      raw_standard_errors[
        number_longitudinal_parameters + 2
      ]
    
    standard_errors[
      number_longitudinal_parameters + 3
    ] <-
      raw_standard_errors[
        number_longitudinal_parameters + 3
      ]
    
    standard_errors[
      number_longitudinal_parameters + 4
    ] <-
      raw_standard_errors[
        number_longitudinal_parameters + 4
      ]
    
    standard_errors[
      number_longitudinal_parameters + 5
    ] <-
      raw_standard_errors[
        number_longitudinal_parameters + 5
      ]
    
    standard_errors[
      number_longitudinal_parameters + 6
    ] <-
      object$Tau *
      raw_standard_errors[
        number_longitudinal_parameters + 6
      ]
  }
  
  z_value <- estimates / standard_errors
  p_value <- 2 * pnorm(
    abs(z_value),
    lower.tail = FALSE
  )
  
  coefficient_table <- data.frame(
    Estimate = estimates,
    Std.Error = standard_errors,
    z.value = z_value,
    p.value = p_value,
    check.names = FALSE
  )
  
  result <- list(
    call = object$call,
    coefficients = coefficient_table,
    logLik = object$loglik,
    objective = object$objective,
    convergence = object$convergence,
    convergence_message =
      object$convergence_message,
    number_subjects =
      nrow(object$survival_data),
    number_longitudinal_records =
      nrow(object$longitudinal_data),
    components_between =
      object$lfpca$KB,
    components_within =
      object$lfpca$KU
  )
  
  class(result) <- "summary.lfjm"
  result
}
