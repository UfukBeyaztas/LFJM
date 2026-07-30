\name{residuals.lfjm}
\alias{residuals.lfjm}
\title{Extract Residuals from an LFJM Fit}
\description{
Extracts conditional response residuals or marginal longitudinal residuals
from a fitted longitudinal functional joint model.
}
\usage{
\method{residuals}{lfjm}(
  object,
  type = c("response", "marginal"),
  ...
)
}
\arguments{
  \item{object}{
  An object of class \code{"lfjm"}, typically returned by
  \code{\link{lfjm}}.
  }

  \item{type}{
  A character string specifying the residual definition. Use
  \code{"response"} for residuals conditional on the fitted
  subject-specific random intercepts, or \code{"marginal"} for residuals
  calculated without those random intercepts. Partial matching is allowed.
  }

  \item{\dots}{
  Additional arguments. They are currently not used.
  }
}
\details{
Let the fitted fixed part of the longitudinal submodel be
\deqn{
\widehat{\mu}_{ij}^{F}
=
\widehat{\beta}_0
+ \widehat{\beta}_{\mathrm{time}}t_{ij}
+ \widehat{\beta}_{Q}Q_i
+ \widehat{\boldsymbol{\xi}}_i^\mathsf{T}
  \widehat{\boldsymbol{\theta}}_B
+ \widehat{\boldsymbol{\zeta}}_{ij}^\mathsf{T}
  \widehat{\boldsymbol{\theta}}_U,
}{
mu_F_hat_ij = beta_0_hat + beta_time_hat t_ij + beta_Q_hat Q_i
+ xi_hat_i' theta_B_hat + zeta_hat_ij' theta_U_hat.
}
The between-subject scores
\eqn{\widehat{\boldsymbol{\xi}}_i} and visit-specific within-subject scores
\eqn{\widehat{\boldsymbol{\zeta}}_{ij}} are obtained from the Stage-1
longitudinal functional principal component analysis and are treated as fixed
quantities during Stage 2.

For \code{type = "response"}, the method returns the conditional response
residual
\deqn{
r_{ij}^{C}
=
Y_{ij}
- \left(\widehat{\mu}_{ij}^{F}+\widehat{u}_i\right),
}{
r_C_ij = Y_ij - (mu_F_hat_ij + u_i_hat),
}
where \eqn{\widehat{u}_i} is the profiled subject-specific random-intercept
mode from the adjusted profile h-likelihood fit. These residuals equal the
observed longitudinal responses minus \code{\link{fitted.lfjm}} values.

For \code{type = "marginal"}, the method omits the fitted random intercept and
returns
\deqn{
r_{ij}^{M}
=
Y_{ij}-\widehat{\mu}_{ij}^{F}
=
r_{ij}^{C}+\widehat{u}_i.
}{
r_M_ij = Y_ij - mu_F_hat_ij = r_C_ij + u_i_hat.
}
The term \dQuote{marginal} therefore refers specifically to exclusion of the
subject-specific random intercept from the longitudinal fitted value. The
functional score contributions remain in
\eqn{\widehat{\mu}_{ij}^{F}}. These are longitudinal residuals; the method
does not calculate martingale, deviance, or Cox--Snell residuals for the
survival submodel.

The residuals are returned in the internally ordered longitudinal-data order
stored in \code{object$longitudinal_data}. During model fitting, this data
frame and the functional-data rows are ordered jointly by subject identifier
and visit time.
}
\value{
A numeric vector with one residual for each row of
\code{object$longitudinal_data}.
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

response_residuals <- residuals(fit)
marginal_residuals <- residuals(
  fit,
  type = "marginal"
)

summary(response_residuals)
summary(marginal_residuals)

random_intercepts <- fit$uhat[
  match(
    fit$longitudinal_data$ID,
    fit$survival_data$ID
  )
]

all.equal(
  marginal_residuals,
  response_residuals + random_intercepts
)
}
}
