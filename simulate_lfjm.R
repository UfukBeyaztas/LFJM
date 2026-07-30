simulate_lfjm <- function(
    setting = c("longitudinal", "baseline"),
    number_subjects = 400,
    number_visits = 8,
    grid = seq(0, 1, length.out = 101),
    visit_grid = seq(0, 5, length.out = number_visits),
    seed = NULL
) {
  setting <- match.arg(setting)
  
  if (setting == "baseline") {
    simulated_data <- simulate_lfjm_core(
      number_subjects = number_subjects,
      number_visits = number_visits,
      grid = grid,
      visit_grid = visit_grid,
      lambda_between = c(3.0, 1.5, 0.6),
      lambda_within = c(1e-8, 1e-8),
      theta_between = c(0.90, -0.55, 0.00),
      theta_within = c(0, 0),
      sigma_u = 0.40,
      sigma_epsilon = 0.30,
      Tau = 1.20,
      gamma0 = -4.20,
      gamma_Q = -0.45,
      alpha = 0.90,
      censoring_rate = 0.05,
      maximum_follow_up = 8.0,
      functional_error_sd = 1e-4,
      slope_scale = 0.0,
      seed = seed
    )
    
    simulated_data$case_type <- "baseline"
    
    functional_columns <-
      functional_column_names("func.X.", simulated_data$grid)
    
    longitudinal_data <- simulated_data$data.long[
      order(
        simulated_data$data.long$ID,
        simulated_data$data.long$obstime
      ),
      ,
      drop = FALSE
    ]
    
    first_visit <- longitudinal_data[
      !duplicated(longitudinal_data$ID),
      c("ID", functional_columns),
      drop = FALSE
    ]
    
    repeated_baseline_curves <- as.matrix(
      first_visit[
        match(longitudinal_data$ID, first_visit$ID),
        functional_columns,
        drop = FALSE
      ]
    )
    
    longitudinal_data[, functional_columns] <- repeated_baseline_curves
    simulated_data$data.long <- longitudinal_data
    simulated_data$true_U_visit[,] <- 0
    simulated_data$truth$theta_U[] <- 0
    simulated_data$truth$beta_U_true[] <- 0
  } else {
    simulated_data <- simulate_lfjm_core(
      number_subjects = number_subjects,
      number_visits = number_visits,
      grid = grid,
      visit_grid = visit_grid,
      lambda_between = c(3.0, 1.5, 0.6),
      lambda_within = c(0.5, 0.25),
      theta_between = c(0.90, -0.55, 0.00),
      theta_within = c(0.45, -0.30),
      sigma_u = 0.40,
      sigma_epsilon = 0.30,
      Tau = 1.20,
      gamma0 = -4.20,
      gamma_Q = -0.45,
      alpha = 0.90,
      censoring_rate = 0.05,
      maximum_follow_up = 8.0,
      functional_error_sd = 0.015,
      slope_scale = 0.0,
      seed = seed
    )
    
    simulated_data$case_type <- "longitudinal"
  }
  
  class(simulated_data) <- c("lfjm_simulation", "list")
  simulated_data
}