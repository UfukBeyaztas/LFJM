coef.lfjm <- function(
    object,
    type = c(
      "all",
      "scalar",
      "between",
      "within",
      "survival",
      "scores"
    ),
    ...
) {
  type <- match.arg(type)
  
  scalar_parameters <- c(
    beta0 = unname(
      object$beta["(Intercept)"]
    ),
    beta_time = unname(
      object$beta["obstime"]
    ),
    beta_Q = unname(
      object$beta["Q1"]
    ),
    sigma_e = object$sigma_e,
    sigma_u = object$sigma_u,
    gamma0 = object$gamma0,
    gamma_Q = object$gamma_Q,
    alpha = object$alpha,
    Tau = object$Tau
  )
  
  if (type == "scalar") {
    return(scalar_parameters)
  }
  
  if (type == "between") {
    return(object$beta_B)
  }
  
  if (type == "within") {
    return(object$beta_U)
  }
  
  if (type == "survival") {
    return(object$surv_B_current)
  }
  
  score_coefficients <- list(
    between = object$theta_B,
    within = object$theta_U
  )
  
  if (type == "scores") {
    return(score_coefficients)
  }
  
  list(
    scalar = scalar_parameters,
    score_coefficients =
      score_coefficients,
    beta_between = object$beta_B,
    beta_within = object$beta_U,
    survival_between =
      object$surv_B_current
  )
}
