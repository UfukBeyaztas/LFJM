\name{fitted.lfjm}
\alias{fitted.lfjm}
\title{Extract Longitudinal Fitted Values from an LFJM Fit}
\description{
Returns the subject-specific fitted values for the longitudinal outcome from a
fitted longitudinal functional joint model.
}
\usage{
\method{fitted}{lfjm}(object, ...)
}
\arguments{
  \item{object}{
    An object of class \code{"lfjm"}, typically returned by
    \code{\link{lfjm}}.
  }
  \item{\dots}{
    Additional arguments. They are currently not used.
  }
}
\details{
The longitudinal functional joint model separates the centered functional
predictor into a persistent subject-level component and a visit-specific
component. For subject \eqn{i} at visit \eqn{j}, the fitted value returned by
this method is
\deqn{
\widehat{\mu}_{ij}
= \widehat{\beta}_0
+ \widehat{\beta}_{\mathrm{time}}t_{ij}
+ \widehat{\beta}_Q Q_i
+ \widehat{\boldsymbol{\xi}}_i^\mathsf{T}
  \widehat{\boldsymbol{\theta}}_B
+ \widehat{\boldsymbol{\zeta}}_{ij}^\mathsf{T}
  \widehat{\boldsymbol{\theta}}_U
+ \widehat{u}_i.
}{
mu_hat_ij = beta0_hat + beta_time_hat t_ij + beta_Q_hat Q_i
+ xi_hat_i' theta_B_hat + zeta_hat_ij' theta_U_hat + u_hat_i.
}
Here, \eqn{\widehat{\boldsymbol{\xi}}_i} and
\eqn{\widehat{\boldsymbol{\zeta}}_{ij}} are the between- and within-level
scores obtained from Stage 1, while
\eqn{\widehat{\boldsymbol{\theta}}_B} and
\eqn{\widehat{\boldsymbol{\theta}}_U} are their Stage-2 coefficient vectors.
The term \eqn{\widehat{u}_i} is the subject-specific random-intercept mode
profiled in the adjusted profile h-likelihood.

These are conditional, subject-specific fitted values. They include both
functional score contributions and the profiled random intercept; therefore,
they differ from marginal fitted values formed from the fixed part alone.
Stage-1 quantities are treated as fixed in accordance with the two-stage
estimation procedure.
}
\value{
A numeric vector containing one subject-specific fitted longitudinal value for
each row of the longitudinal data stored in \code{object}. The values follow
the row ordering used in the fitted model, namely subject identifier followed
by visit time.
}
\examples{
\dontrun{
simulated_data <- simulate_lfjm(
  setting = "longitudinal",
  number_subjects = 100,
  number_visits = 6
)

functional_columns <- paste0(
  "func.X.",
  seq_along(simulated_data$grid)
)

functional_data <- as.matrix(
  simulated_data$data.long[, functional_columns, drop = FALSE]
)

fit <- lfjm(
  longitudinal_data = simulated_data$data.long,
  survival_data = simulated_data$data.surv,
  functional_data = functional_data,
  grid = simulated_data$grid
)

fitted_values <- fitted(fit)
head(fitted_values)
}
}
