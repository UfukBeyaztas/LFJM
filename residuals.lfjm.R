residuals.lfjm <- function(
    object,
    type = c("response", "marginal"),
    ...
) {
  type <- match.arg(type)
  
  if (type == "response") {
    return(
      object$longitudinal_data$Y -
        fitted(object)
    )
  }
  
  object$longitudinal_data$Y -
    longitudinal_fixed_part(object)
}