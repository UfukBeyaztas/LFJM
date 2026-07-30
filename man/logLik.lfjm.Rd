\name{logLik.lfjm}
\alias{logLik.lfjm}
\title{Extract the Adjusted Log-Likelihood from an LFJM Fit}
\description{
Extracts the maximized adjusted profile h-likelihood criterion from a fitted
longitudinal functional joint model and returns it as an object of class
\code{"logLik"}.
}
\usage{
\method{logLik}{lfjm}(object, ...)
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
The longitudinal functional joint model contains a subject-specific random
intercept that is shared by the longitudinal and survival submodels. Let
\eqn{h_i(u_i;\boldsymbol{\psi})} denote the hierarchical log-likelihood
contribution of subject \eqn{i}, including its longitudinal, survival, and
random-effect terms. For fixed model parameters
\eqn{\boldsymbol{\psi}}, the random-intercept mode is
\deqn{
\widehat{u}_i =
\mathop{\mathrm{arg\,max}}_{u_i}
h_i(u_i;\boldsymbol{\psi}).
}{
u_hat_i = arg max over u_i of h_i(u_i; psi).
}

The fitted criterion is the adjusted profile h-likelihood
\deqn{
\widetilde{h}(\boldsymbol{\psi})
=
\sum_{i=1}^{N}
\left[
h_i(\widehat{u}_i;\boldsymbol{\psi})
-\frac{1}{2}
\log\left\{\frac{H_i}{2\pi}\right\}
\right],
}{
h_tilde(psi) =
sum_i [h_i(u_hat_i; psi) - 0.5 log\{H_i/(2 pi)\}].
}
Here,
\eqn{H_i=-\partial^2h_i(u_i;\boldsymbol{\psi})/\partial u_i^2} evaluated at
\eqn{\widehat{u}_i} is the positive curvature of the hierarchical
log-likelihood. In the implemented random-intercept model,
\deqn{
H_i =
\frac{n_i}{\sigma_e^2}
+ \alpha^2\Lambda_i
+ \frac{1}{\sigma_u^2},
}{
H_i = n_i / sigma_e^2 + alpha^2 Lambda_i + 1 / sigma_u^2,
}
where \eqn{n_i} is the number of longitudinal observations and
\eqn{\Lambda_i} is the fitted cumulative hazard contribution for subject
\eqn{i}.

This adjustment is the first-order Laplace approximation to the marginal
log-likelihood obtained by integrating out the random intercept:
\deqn{
\log\int
\exp\{h_i(u_i;\boldsymbol{\psi})\}\,du_i
=
h_i(\widehat{u}_i;\boldsymbol{\psi})
-\frac{1}{2}\log\left\{\frac{H_i}{2\pi}\right\}
+R_i.
}{
log integral exp\{h_i(u_i; psi)\} du_i =
h_i(u_hat_i; psi) - 0.5 log\{H_i/(2 pi)\} + R_i.
}
The returned criterion omits the higher-order remainder \eqn{R_i}. When the
association parameter \eqn{\alpha=0}, the hierarchical log-likelihood is
quadratic in the random intercept, \eqn{R_i=0}, and the adjustment is exact.
For nonzero association, the returned value is an adjusted profile
h-likelihood approximation rather than an exactly marginalized
log-likelihood.

The \code{"df"} attribute equals \code{object$npar}, the length of the fitted
Stage-2 parameter vector. It includes the longitudinal score coefficients and
finite-dimensional model parameters but excludes the profiled
subject-specific random intercepts and all Stage-1 LFPCA quantities. The
\code{"nobs"} attribute equals the number of rows in
\code{object$longitudinal_data}; it therefore counts longitudinal visit-level
observations rather than independent subjects or survival records.
}
\value{
A numeric value of class \code{"logLik"} containing the maximized adjusted
profile h-likelihood stored in \code{object$loglik}. It has the following
attributes:
\describe{
  \item{\code{df}}{
  The number of fitted Stage-2 parameters, obtained from \code{object$npar}.
  }
  \item{\code{nobs}}{
  The number of rows in the longitudinal data stored in the fitted object.
  }
}
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

adjusted_loglikelihood <- logLik(fit)
adjusted_loglikelihood
as.numeric(adjusted_loglikelihood)
attr(adjusted_loglikelihood, "df")
attr(adjusted_loglikelihood, "nobs")

AIC(fit)
BIC(fit)
}
}
