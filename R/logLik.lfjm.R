logLik.lfjm <- function(object, ...) {
  structure(
    object$loglik,
    df = object$npar,
    nobs = nrow(object$longitudinal_data),
    class = "logLik"
  )
}
