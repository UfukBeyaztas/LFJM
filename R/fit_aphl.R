fit_aphl <- function(
    longitudinal_data,
    survival_data,
    lfpca,
    grid,
    quadrature_points = 32,
    max_iterations = 500,
    multistart = TRUE,
    standard_errors = TRUE
) {
  longitudinal_data <- longitudinal_data[
    order(
      longitudinal_data$ID,
      longitudinal_data$obstime
    ),
    ,
    drop = FALSE
  ]
  
  subject_ids <- sort(unique(longitudinal_data$ID))
  between_scores <- lfpca$scoreB_subject[
    as.character(subject_ids),
    ,
    drop = FALSE
  ]
  
  between_scores_long <- between_scores[
    match(longitudinal_data$ID, subject_ids),
    ,
    drop = FALSE
  ]
  
  within_scores_long <- lfpca$scoreU_visit
  number_between_components <- ncol(between_scores)
  number_within_components <- ncol(within_scores_long)
  
  between_names <- if (number_between_components > 0) {
    paste0("Bscore", seq_len(number_between_components))
  } else {
    character(0)
  }
  
  within_names <- if (number_within_components > 0) {
    paste0("Uscore", seq_len(number_within_components))
  } else {
    character(0)
  }
  
  colnames(between_scores_long) <- between_names
  colnames(within_scores_long) <- within_names
  
  longitudinal_design <- as.matrix(
    cbind(
      "(Intercept)" = 1,
      obstime = longitudinal_data$obstime,
      Q1 = longitudinal_data$Q1,
      between_scores_long,
      within_scores_long
    )
  )
  
  response <- as.numeric(longitudinal_data$Y)
  subject_index <- match(longitudinal_data$ID, subject_ids)
  number_subjects <- length(subject_ids)
  number_longitudinal_parameters <- ncol(longitudinal_design)
  
  parameter_index <- list(
    intercept = 1,
    obstime = 2,
    Q1 = 3,
    theta_between =
      3 + seq_len(number_between_components),
    theta_within =
      3 +
      number_between_components +
      seq_len(number_within_components)
  )
  
  observations_per_subject <- as.numeric(
    tabulate(
      subject_index,
      nbins = number_subjects
    )
  )
  
  survival_data <- survival_data[
    match(subject_ids, survival_data$ID),
    ,
    drop = FALSE
  ]
  
  survival_covariate <- survival_data$Q1
  survival_time <- pmax(survival_data$time, 1e-8)
  event <- as.numeric(survival_data$event)
  quadrature <- gauss_legendre(quadrature_points)
  
  evaluate_criterion <- function(parameters) {
    beta <- parameters[
      seq_len(number_longitudinal_parameters)
    ]
    residual_variance <-
      exp(
        2 *
          parameters[number_longitudinal_parameters + 1]
      )
    random_intercept_variance <-
      exp(
        2 *
          parameters[number_longitudinal_parameters + 2]
      )
    gamma0 <-
      parameters[number_longitudinal_parameters + 3]
    gamma_Q <-
      parameters[number_longitudinal_parameters + 4]
    alpha <-
      parameters[number_longitudinal_parameters + 5]
    Tau <-
      exp(
        parameters[number_longitudinal_parameters + 6]
      )
    
    residual <- as.numeric(
      response -
        longitudinal_design %*% beta
    )
    
    residual_sum <- rowsum(
      residual,
      subject_index
    )[, 1]
    squared_residual_sum <- rowsum(
      residual^2,
      subject_index
    )[, 1]
    
    theta_between <- beta[
      parameter_index$theta_between
    ]
    
    current_value_intercept <-
      beta[1] +
      beta[parameter_index$Q1] *
      survival_covariate +
      as.numeric(
        between_scores %*%
          theta_between
      )
    
    current_value_slope <-
      beta[parameter_index$obstime]
    
    survival_linear_predictor <-
      gamma0 +
      gamma_Q * survival_covariate +
      alpha * current_value_intercept
    
    transformed_integral <-
      survival_time^Tau *
      as.numeric(
        bounded_exponential(
          outer(
            alpha *
              current_value_slope *
              survival_time,
            quadrature$nodes^(1 / Tau)
          )
        ) %*%
          quadrature$weights
      )
    
    cumulative_hazard_without_random_effect <-
      bounded_exponential(
        survival_linear_predictor
      ) *
      transformed_integral
    
    random_intercept <- rep(0, number_subjects)
    
    for (iteration in seq_len(60)) {
      cumulative_hazard <-
        bounded_exponential(
          alpha * random_intercept
        ) *
        cumulative_hazard_without_random_effect
      
      score <-
        (
          residual_sum -
            observations_per_subject *
            random_intercept
        ) /
        residual_variance +
        event * alpha -
        alpha * cumulative_hazard -
        random_intercept /
        random_intercept_variance
      
      second_derivative <-
        -observations_per_subject /
        residual_variance -
        alpha^2 * cumulative_hazard -
        1 / random_intercept_variance
      
      step <- -score / second_derivative
      step[!is.finite(step)] <- 0
      step <- pmax(pmin(step, 2), -2)
      
      random_intercept <-
        pmax(
          pmin(
            random_intercept + step,
            10
          ),
          -10
        )
      
      if (max(abs(step)) < 1e-9) {
        break
      }
    }
    
    cumulative_hazard <-
      bounded_exponential(
        alpha * random_intercept
      ) *
      cumulative_hazard_without_random_effect
    
    second_derivative <-
      -observations_per_subject /
      residual_variance -
      alpha^2 * cumulative_hazard -
      1 / random_intercept_variance
    
    residual_sum_of_squares <-
      squared_residual_sum -
      2 * random_intercept * residual_sum +
      observations_per_subject *
      random_intercept^2
    
    longitudinal_loglikelihood <-
      -observations_per_subject /
      2 *
      log(2 * pi * residual_variance) -
      residual_sum_of_squares /
      (2 * residual_variance)
    
    log_hazard <-
      log(Tau) +
      (Tau - 1) *
      log(survival_time) +
      survival_linear_predictor +
      alpha * random_intercept +
      alpha *
      current_value_slope *
      survival_time
    
    survival_loglikelihood <-
      event * log_hazard -
      cumulative_hazard
    
    random_effect_loglikelihood <-
      -0.5 *
      log(
        2 *
          pi *
          random_intercept_variance
      ) -
      random_intercept^2 /
      (
        2 *
          random_intercept_variance
      )
    
    list(
      marginal =
        longitudinal_loglikelihood +
        survival_loglikelihood +
        random_effect_loglikelihood +
        0.5 * log(2 * pi) -
        0.5 * log(-second_derivative),
      random_intercept = random_intercept,
      risk =
        survival_linear_predictor +
        alpha * random_intercept +
        alpha *
        current_value_slope *
        survival_time
    )
  }
  
  objective <- function(parameters) {
    value <- -sum(
      evaluate_criterion(parameters)$marginal
    )
    
    if (!is.finite(value)) {
      return(1e12)
    }
    
    value
  }
  
  initial_beta <- as.numeric(
    solve_system(
      crossprod(longitudinal_design) +
        1e-8 *
        diag(number_longitudinal_parameters),
      crossprod(longitudinal_design, response)
    )
  )
  initial_beta[!is.finite(initial_beta)] <- 0
  
  initial_residual <- as.numeric(
    response -
      longitudinal_design %*%
      initial_beta
  )
  
  initial_subject_mean <- as.numeric(
    tapply(
      initial_residual,
      subject_index,
      mean
    )
  )
  
  initial_residual_sd <- sqrt(
    max(
      mean(
        (
          initial_residual -
            initial_subject_mean[subject_index]
        )^2
      ),
      1e-3
    )
  )
  
  initial_survival_fit <- tryCatch(
    survreg(
      Surv(
        survival_time,
        event
      ) ~
        survival_covariate,
      dist = "weibull"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(initial_survival_fit)) {
    initial_scale <- initial_survival_fit$scale
    initial_survival_coefficients <-
      coef(initial_survival_fit)
    
    initial_log_Tau <- log(
      clamp_values(
        1 / initial_scale,
        0.3,
        4
      )
    )
    
    initial_gamma0 <-
      -initial_survival_coefficients[1] /
      initial_scale
    
    initial_gamma_Q <- if (
      "survival_covariate" %in%
      names(initial_survival_coefficients)
    ) {
      -initial_survival_coefficients[
        "survival_covariate"
      ] /
        initial_scale
    } else {
      0
    }
  } else {
    initial_log_Tau <- 0
    initial_gamma0 <- -3
    initial_gamma_Q <- 0
  }
  
  initial_parameters <- c(
    initial_beta,
    log(initial_residual_sd),
    log(
      max(
        sd(initial_subject_mean),
        0.1
      )
    ),
    initial_gamma0,
    initial_gamma_Q,
    0.1,
    initial_log_Tau
  )
  initial_parameters[!is.finite(initial_parameters)] <- 0
  
  lower_bounds <- c(
    rep(
      -50,
      number_longitudinal_parameters
    ),
    log(1e-2),
    log(1e-2),
    -12,
    -8,
    -3,
    log(0.2)
  )
  
  upper_bounds <- c(
    rep(
      50,
      number_longitudinal_parameters
    ),
    log(20),
    log(20),
    12,
    8,
    3,
    log(5)
  )
  
  initial_parameters <- pmin(
    pmax(
      initial_parameters,
      lower_bounds + 1e-6
    ),
    upper_bounds - 1e-6
  )
  
  optimize_aphl <- function(starting_values) {
    optim(
      starting_values,
      objective,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(
        maxit = max_iterations,
        factr = 1e7
      )
    )
  }
  
  alpha_starting_values <- if (multistart) {
    c(-0.5, 0, 0.5, 1)
  } else {
    0.1
  }
  
  candidate_fits <- lapply(
    alpha_starting_values,
    function(initial_alpha) {
      candidate_parameters <- initial_parameters
      candidate_parameters[
        number_longitudinal_parameters + 5
      ] <- initial_alpha
      
      tryCatch(
        optimize_aphl(candidate_parameters),
        error = function(e) NULL
      )
    }
  )
  
  candidate_fits <- candidate_fits[
    !vapply(
      candidate_fits,
      is.null,
      logical(1)
    )
  ]
  
  if (length(candidate_fits) == 0) {
    stop("All APHL optimization starts failed.")
  }
  
  fitted_model <- candidate_fits[[
    which.min(
      vapply(
        candidate_fits,
        function(candidate) {
          candidate$value
        },
        numeric(1)
      )
    )
  ]]
  
  if (multistart) {
    Nelder_Mead_fit <- tryCatch(
      optim(
        fitted_model$par,
        objective,
        method = "Nelder-Mead",
        control = list(
          maxit = 2000,
          reltol = 1e-9
        )
      ),
      error = function(e) NULL
    )
    
    if (!is.null(Nelder_Mead_fit)) {
      polished_start <- pmin(
        pmax(
          Nelder_Mead_fit$par,
          lower_bounds + 1e-6
        ),
        upper_bounds - 1e-6
      )
      
      polished_fit <- tryCatch(
        optimize_aphl(polished_start),
        error = function(e) NULL
      )
      
      if (!is.null(polished_fit) &&
          polished_fit$value < fitted_model$value) {
        fitted_model <- polished_fit
      }
    }
  }
  
  fitted_parameters <- fitted_model$par
  beta <- fitted_parameters[
    seq_len(number_longitudinal_parameters)
  ]
  names(beta) <- colnames(longitudinal_design)
  
  fitted_criterion <- evaluate_criterion(
    fitted_parameters
  )
  
  between_map <- dual_basis_map(
    lfpca$phi0,
    grid
  )
  within_map <- dual_basis_map(
    lfpca$phiU,
    grid
  )
  
  beta_between <- as.numeric(
    between_map %*%
      beta[parameter_index$theta_between]
  )
  beta_within <- as.numeric(
    within_map %*%
      beta[parameter_index$theta_within]
  )
  alpha <- fitted_parameters[
    number_longitudinal_parameters + 5
  ]
  
  parameter_names <- c(
    colnames(longitudinal_design),
    "log_sigma_e",
    "log_sigma_u",
    "gamma0",
    "gamma_Q",
    "alpha",
    "log_Tau"
  )
  
  names(fitted_parameters) <- parameter_names
  
  covariance_matrix <- NULL
  se_beta_between <- NULL
  se_beta_within <- NULL
  se_survival_between <- NULL
  
  if (standard_errors) {
    Hessian <- tryCatch(
      optimHess(
        fitted_parameters,
        objective
      ),
      error = function(e) NULL
    )
    
    covariance_matrix <- if (!is.null(Hessian)) {
      tryCatch(
        solve(Hessian),
        error = function(e) {
          tryCatch(
            ginv(Hessian),
            error = function(e) NULL
          )
        }
      )
    } else {
      NULL
    }
    
    if (!is.null(covariance_matrix)) {
      dimnames(covariance_matrix) <- list(
        parameter_names,
        parameter_names
      )
      
      between_index <-
        parameter_index$theta_between
      within_index <-
        parameter_index$theta_within
      alpha_index <-
        number_longitudinal_parameters + 5
      
      se_beta_between <- sqrt(
        pmax(
          rowSums(
            (
              between_map %*%
                covariance_matrix[
                  between_index,
                  between_index,
                  drop = FALSE
                ]
            ) *
              between_map
          ),
          0
        )
      )
      
      se_beta_within <- sqrt(
        pmax(
          rowSums(
            (
              within_map %*%
                covariance_matrix[
                  within_index,
                  within_index,
                  drop = FALSE
                ]
            ) *
              within_map
          ),
          0
        )
      )
      
      survival_derivative <- cbind(
        beta_between,
        alpha * between_map
      )
      
      covariance_alpha_between <-
        covariance_matrix[
          c(alpha_index, between_index),
          c(alpha_index, between_index),
          drop = FALSE
        ]
      
      se_survival_between <- sqrt(
        pmax(
          rowSums(
            (
              survival_derivative %*%
                covariance_alpha_between
            ) *
              survival_derivative
          ),
          0
        )
      )
    }
  }
  
  list(
    method = "Proposed LFPCA-APHL",
    beta = beta,
    theta_B =
      beta[parameter_index$theta_between],
    theta_U =
      beta[parameter_index$theta_within],
    beta_B = beta_between,
    beta_U = beta_within,
    surv_B_current = alpha * beta_between,
    se_beta_B = se_beta_between,
    se_beta_U = se_beta_within,
    se_surv_B = se_survival_between,
    sigma_e = exp(
      fitted_parameters[
        number_longitudinal_parameters + 1
      ]
    ),
    sigma_u = exp(
      fitted_parameters[
        number_longitudinal_parameters + 2
      ]
    ),
    gamma0 = fitted_parameters[
      number_longitudinal_parameters + 3
    ],
    gamma_Q = fitted_parameters[
      number_longitudinal_parameters + 4
    ],
    alpha = alpha,
    Tau = exp(
      fitted_parameters[
        number_longitudinal_parameters + 6
      ]
    ),
    convergence = fitted_model$convergence,
    convergence_message = fitted_model$message,
    uhat = fitted_criterion$random_intercept,
    risk = fitted_criterion$risk,
    loglik = sum(fitted_criterion$marginal),
    npar = length(fitted_parameters),
    parameter_vector = fitted_parameters,
    covariance_matrix = covariance_matrix,
    objective = fitted_model$value,
    lfpca = lfpca
  )
}
