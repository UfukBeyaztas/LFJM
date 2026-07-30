predict.lfjm <- function(
    object,
    type = c(
      "longitudinal",
      "current_value",
      "risk",
      "random_effects",
      "scores"
    ),
    landmark = NULL,
    ...
) {
  type <- match.arg(type)
  
  if (type == "longitudinal") {
    fixed_part <- longitudinal_fixed_part(object)
    
    fitted_values <-
      fixed_part +
      object$uhat[
        match(
          object$longitudinal_data$ID,
          object$survival_data$ID
        )
      ]
    
    return(as.numeric(fitted_values))
  }
  
  if (type == "scores") {
    return(
      list(
        between = object$lfpca$scoreB_subject,
        within = object$lfpca$scoreU_visit
      )
    )
  }
  
  if (type == "random_effects" &&
      is.null(landmark)) {
    predicted_random_intercepts <- object$uhat
    names(predicted_random_intercepts) <-
      object$survival_data$ID
    return(predicted_random_intercepts)
  }
  
  if (is.null(landmark) ||
      !is.numeric(landmark) ||
      length(landmark) != 1 ||
      !is.finite(landmark)) {
    stop(
      "'landmark' must be supplied for current-value or risk prediction."
    )
  }
  
  predictions <- landmark_prediction(
    object,
    landmark
  )
  
  if (type == "current_value") {
    values <- predictions$current_value
    names(values) <- predictions$ID
    return(values)
  }
  
  if (type == "risk") {
    values <- predictions$risk
    names(values) <- predictions$ID
    return(values)
  }
  
  values <- predictions$random_intercept
  names(values) <- predictions$ID
  values
}