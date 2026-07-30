vcov.lfjm <- function(object, ...) {
  if (is.null(object$covariance_matrix)) {
    stop(
      "The covariance matrix was not computed for this fit."
    )
  }
  
  object$covariance_matrix
}