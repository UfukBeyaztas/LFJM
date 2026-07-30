simulate_lfjm_core <- function(
    number_subjects = 400,
    number_visits = 8,
    grid = seq(0, 1, length.out = 101),
    visit_grid = seq(0, 5, length.out = number_visits),
    beta0 = 1.0,
    beta_time = 0.28,
    beta_Q = 0.45,
    sigma_u = 0.40,
    sigma_epsilon = 0.30,
    Tau = 1.20,
    gamma0 = -4.20,
    gamma_Q = -0.45,
    alpha = 0.90,
    censoring_rate = 0.05,
    maximum_follow_up = 8.0,
    lambda_between = c(3.0, 1.5, 0.6),
    lambda_within = c(0.5, 0.25),
    theta_between = c(0.90, -0.55, 0.00),
    theta_within = c(0.45, -0.30),
    functional_error_sd = 0.015,
    slope_scale = 0.0,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (length(visit_grid) != number_visits) {
    stop("'visit_grid' must have length equal to 'number_visits'.")
  }
  
  grid <- as.numeric(grid)
  number_grid_points <- length(grid)
  number_between_components <- length(lambda_between)
  number_within_components <- length(lambda_within)
  
  theta_between <- theta_between[seq_len(number_between_components)]
  theta_within <- theta_within[seq_len(number_within_components)]
  
  time_mean <- mean(visit_grid)
  time_sd <- sd(rep(visit_grid, number_subjects))
  
  if (!is.finite(time_sd) || time_sd < 1e-8) {
    time_sd <- 1
  }
  
  standardize_time <- function(time) {
    (time - time_mean) / time_sd
  }
  
  phi0 <- make_latent_basis(
    grid,
    number_between_components,
    "B0"
  )
  phi1 <- make_latent_basis(
    grid,
    number_between_components,
    "B1"
  )
  phiU <- make_latent_basis(
    grid,
    number_within_components,
    "U"
  )
  
  beta_between_true <- as.numeric(phi0 %*% theta_between)
  beta_within_true <- as.numeric(phiU %*% theta_within)
  
  mean_function <- function(standardized_time) {
    0.25 * sin(pi * grid) +
      0.10 * standardized_time * cos(2 * pi * grid)
  }
  
  Q <- rbinom(number_subjects, 1, 0.5)
  random_intercept <- rnorm(number_subjects, 0, sigma_u)
  
  between_scores <- matrix(
    rnorm(number_subjects * number_between_components),
    number_subjects,
    number_between_components
  ) %*%
    diag(
      sqrt(lambda_between),
      number_between_components,
      number_between_components
    )
  
  between_intercept <- between_scores %*% t(phi0)
  between_slope <- slope_scale * (between_scores %*% t(phi1))
  between_effect <- as.numeric(between_scores %*% theta_between)
  
  current_value <- function(subject, time) {
    time <- as.numeric(time)
    beta0 +
      beta_time * time +
      beta_Q * Q[subject] +
      between_effect[subject] +
      random_intercept[subject]
  }
  
  event_time <- numeric(number_subjects)
  
  for (subject in seq_len(number_subjects)) {
    exponential_draw <- -log(runif(1))
    
    hazard <- function(time) {
      time <- pmax(time, 1e-8)
      linear_predictor <-
        gamma0 +
        gamma_Q * Q[subject] +
        alpha * current_value(subject, time)
      
      Tau *
        time^(Tau - 1) *
        exp(limit_values(linear_predictor, 35))
    }
    
    cumulative_hazard <- function(time) {
      if (time <= 1e-8) {
        return(0)
      }
      
      integrate(
        hazard,
        1e-8,
        time,
        rel.tol = 1e-7,
        subdivisions = 100
      )$value
    }
    
    root_function <- function(time) {
      cumulative_hazard(time) - exponential_draw
    }
    
    upper <- maximum_follow_up
    function_at_upper <- tryCatch(
      root_function(upper),
      error = function(e) NA_real_
    )
    
    while (is.finite(function_at_upper) &&
           function_at_upper < 0 &&
           upper < 60) {
      upper <- 1.5 * upper
      function_at_upper <- tryCatch(
        root_function(upper),
        error = function(e) NA_real_
      )
    }
    
    event_time[subject] <- if (!is.finite(function_at_upper) ||
                               function_at_upper < 0) {
      upper
    } else {
      uniroot(
        root_function,
        c(1e-8, upper),
        tol = 1e-7
      )$root
    }
  }
  
  censoring_time <- pmin(
    rexp(number_subjects, censoring_rate),
    maximum_follow_up
  )
  
  observed_time <- pmin(event_time, censoring_time)
  event <- as.integer(event_time <= censoring_time)
  
  longitudinal_list <- vector("list", number_subjects)
  functional_data <- matrix(
    NA_real_,
    number_subjects * number_visits,
    number_grid_points
  )
  true_between_visit <- matrix(
    NA_real_,
    number_subjects * number_visits,
    number_grid_points
  )
  true_within_visit <- matrix(
    NA_real_,
    number_subjects * number_visits,
    number_grid_points
  )
  
  row_counter <- 0
  
  for (subject in seq_len(number_subjects)) {
    subject_rows <- vector("list", number_visits)
    
    for (visit in seq_len(number_visits)) {
      row_counter <- row_counter + 1
      visit_time <- visit_grid[visit]
      standardized_time <- standardize_time(visit_time)
      
      within_scores <- rnorm(
        number_within_components,
        0,
        sqrt(lambda_within)
      )
      
      within_curve <- as.numeric(within_scores %*% t(phiU))
      between_curve <-
        between_intercept[subject, ] +
        standardized_time * between_slope[subject, ]
      
      functional_data[row_counter, ] <-
        mean_function(standardized_time) +
        between_curve +
        within_curve +
        rnorm(number_grid_points, 0, functional_error_sd)
      
      true_between_visit[row_counter, ] <- between_curve
      true_within_visit[row_counter, ] <- within_curve
      
      conditional_mean <-
        beta0 +
        beta_time * visit_time +
        beta_Q * Q[subject] +
        between_effect[subject] +
        random_intercept[subject] +
        sum(theta_within * within_scores)
      
      subject_rows[[visit]] <- data.frame(
        ID = subject,
        obstime = visit_time,
        row_id = row_counter,
        Q1 = Q[subject],
        time = observed_time[subject],
        event = event[subject],
        Y = conditional_mean + rnorm(1, 0, sigma_epsilon),
        true_mu_BU = conditional_mean
      )
    }
    
    longitudinal_list[[subject]] <- do.call(rbind, subject_rows)
  }
  
  longitudinal_data <- do.call(rbind, longitudinal_list)
  colnames(functional_data) <-
    functional_column_names("func.X.", grid)
  
  longitudinal_data <- cbind(
    longitudinal_data,
    as.data.frame(functional_data)
  )
  
  survival_data <- data.frame(
    ID = seq_len(number_subjects),
    time = observed_time,
    event = event,
    Q1 = Q,
    true_u = random_intercept
  )
  
  follow_up_by_subject <- survival_data$time
  names(follow_up_by_subject) <- as.character(survival_data$ID)
  
  keep <-
    longitudinal_data$obstime <=
    follow_up_by_subject[as.character(longitudinal_data$ID)] + 1e-10
  keep[is.na(keep)] <- FALSE
  
  longitudinal_data <- longitudinal_data[keep, , drop = FALSE]
  longitudinal_data <- longitudinal_data[
    order(longitudinal_data$ID, longitudinal_data$obstime),
    ,
    drop = FALSE
  ]
  survival_data <- survival_data[
    order(survival_data$ID),
    ,
    drop = FALSE
  ]
  
  retained_rows <- longitudinal_data$row_id
  
  list(
    data.long = longitudinal_data,
    data.surv = survival_data,
    grid = grid,
    sgrid = grid,
    true_B_visit = true_between_visit[retained_rows, , drop = FALSE],
    true_U_visit = true_within_visit[retained_rows, , drop = FALSE],
    truth = list(
      beta0 = beta0,
      beta_time = beta_time,
      beta_Q = beta_Q,
      sigma_u = sigma_u,
      sigma_eps = sigma_epsilon,
      Tau = Tau,
      gamma0 = gamma0,
      gamma_Q = gamma_Q,
      alpha = alpha,
      theta_B = theta_between,
      theta_U = theta_within,
      beta_B_true = beta_between_true,
      beta_U_true = beta_within_true,
      surv_B_current_true = alpha * beta_between_true,
      phi0 = phi0,
      phiU = phiU,
      xi = between_scores
    )
  )
}