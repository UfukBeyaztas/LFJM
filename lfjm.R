lfjm <- function(
    longitudinal_data,
    survival_data,
    functional_data,
    grid,
    id = "ID",
    visit_time = "obstime",
    response = "Y",
    longitudinal_covariate = "Q1",
    survival_time = "time",
    event = "event",
    survival_covariate = "Q1",
    control = lfjm_control()
) {
  if (!inherits(control, "lfjm_control")) {
    stop("'control' must be created by lfjm_control().")
  }
  
  validate_lfjm_data(
    longitudinal_data = longitudinal_data,
    survival_data = survival_data,
    functional_data = functional_data,
    grid = grid,
    id = id,
    visit_time = visit_time,
    response = response,
    longitudinal_covariate = longitudinal_covariate,
    survival_time = survival_time,
    event = event,
    survival_covariate = survival_covariate
  )
  
  longitudinal_order <- order(
    longitudinal_data[[id]],
    longitudinal_data[[visit_time]]
  )
  
  longitudinal_data <- longitudinal_data[
    longitudinal_order,
    ,
    drop = FALSE
  ]
  functional_data <- as.matrix(functional_data)[
    longitudinal_order,
    ,
    drop = FALSE
  ]
  
  model_longitudinal_data <- data.frame(
    ID = longitudinal_data[[id]],
    obstime = as.numeric(
      longitudinal_data[[visit_time]]
    ),
    Q1 = as.numeric(
      longitudinal_data[[longitudinal_covariate]]
    ),
    Y = as.numeric(
      longitudinal_data[[response]]
    ),
    stringsAsFactors = FALSE
  )
  
  subject_ids <- sort(
    unique(model_longitudinal_data$ID)
  )
  survival_order <- match(
    subject_ids,
    survival_data[[id]]
  )
  
  survival_data <- survival_data[
    survival_order,
    ,
    drop = FALSE
  ]
  
  model_survival_data <- data.frame(
    ID = survival_data[[id]],
    time = as.numeric(
      survival_data[[survival_time]]
    ),
    event = as.numeric(
      survival_data[[event]]
    ),
    Q1 = as.numeric(
      survival_data[[survival_covariate]]
    ),
    stringsAsFactors = FALSE
  )
  
  lfpca_fit <- lfjm_lfpca(
    functional_data = functional_data,
    subject = model_longitudinal_data$ID,
    time = model_longitudinal_data$obstime,
    grid = grid,
    pve_between = control$pve_between,
    pve_within = control$pve_within,
    components_between =
      control$components_between,
    components_within =
      control$components_within,
    mean_degree = control$mean_degree,
    presmooth = control$presmooth,
    smoothing_df = control$smoothing_df,
    use_functional_slope =
      control$use_functional_slope
  )
  
  fitted_model <- fit_aphl(
    longitudinal_data = model_longitudinal_data,
    survival_data = model_survival_data,
    lfpca = lfpca_fit,
    grid = grid,
    quadrature_points =
      control$quadrature_points,
    max_iterations =
      control$max_iterations,
    multistart =
      control$multistart,
    standard_errors =
      control$standard_errors
  )
  
  fitted_model$call <- match.call()
  fitted_model$longitudinal_data <-
    model_longitudinal_data
  fitted_model$survival_data <-
    model_survival_data
  fitted_model$functional_data <-
    functional_data
  fitted_model$grid <- as.numeric(grid)
  fitted_model$control <- control
  fitted_model$variable_names <- list(
    id = id,
    visit_time = visit_time,
    response = response,
    longitudinal_covariate =
      longitudinal_covariate,
    survival_time = survival_time,
    event = event,
    survival_covariate =
      survival_covariate
  )
  
  class(fitted_model) <- "lfjm"
  fitted_model
}