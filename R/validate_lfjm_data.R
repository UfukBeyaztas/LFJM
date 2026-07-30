validate_lfjm_data <- function(
    longitudinal_data,
    survival_data,
    functional_data,
    grid,
    id,
    visit_time,
    response,
    longitudinal_covariate,
    survival_time,
    event,
    survival_covariate
) {
  required_longitudinal_columns <- c(
    id,
    visit_time,
    response,
    longitudinal_covariate
  )
  
  missing_longitudinal_columns <- setdiff(
    required_longitudinal_columns,
    names(longitudinal_data)
  )
  
  if (length(missing_longitudinal_columns) > 0) {
    stop(
      "Missing longitudinal columns: ",
      paste(
        missing_longitudinal_columns,
        collapse = ", "
      )
    )
  }
  
  required_survival_columns <- c(
    id,
    survival_time,
    event,
    survival_covariate
  )
  
  missing_survival_columns <- setdiff(
    required_survival_columns,
    names(survival_data)
  )
  
  if (length(missing_survival_columns) > 0) {
    stop(
      "Missing survival columns: ",
      paste(
        missing_survival_columns,
        collapse = ", "
      )
    )
  }
  
  functional_data <- as.matrix(functional_data)
  
  if (nrow(functional_data) != nrow(longitudinal_data)) {
    stop(
      "'functional_data' must have one row for each longitudinal record."
    )
  }
  
  if (ncol(functional_data) != length(grid)) {
    stop(
      "'functional_data' must have one column for each functional grid point."
    )
  }
  
  trapezoidal_weights(grid)
  
  if (anyDuplicated(survival_data[[id]])) {
    stop(
      "'survival_data' must contain exactly one row per subject."
    )
  }
  
  longitudinal_subjects <-
    sort(unique(longitudinal_data[[id]]))
  survival_subjects <-
    sort(unique(survival_data[[id]]))
  
  if (!identical(
    as.character(longitudinal_subjects),
    as.character(survival_subjects)
  )) {
    stop(
      "The subject identifiers in the longitudinal and survival data must match."
    )
  }
  
  if (any(!is.finite(longitudinal_data[[visit_time]])) ||
      any(!is.finite(longitudinal_data[[response]])) ||
      any(!is.finite(
        longitudinal_data[[longitudinal_covariate]]
      ))) {
    stop(
      "The longitudinal time, response, and scalar covariate must be finite."
    )
  }
  
  if (any(!is.finite(survival_data[[survival_time]])) ||
      any(survival_data[[survival_time]] <= 0)) {
    stop("All observed survival times must be positive and finite.")
  }
  
  event_values <- survival_data[[event]]
  
  if (any(!is.finite(event_values)) ||
      any(!event_values %in% c(0, 1))) {
    stop("The event indicator must contain only 0 and 1.")
  }
  
  if (any(!is.finite(
    survival_data[[survival_covariate]]
  ))) {
    stop("The survival covariate must be finite.")
  }
  
  if (any(!is.finite(functional_data))) {
    stop("'functional_data' must contain only finite values.")
  }
  
  invisible(TRUE)
}
