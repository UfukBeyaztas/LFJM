\name{vcov.lfjm}
\alias{vcov.lfjm}
\title{Extract the Covariance Matrix from an LFJM Fit}
\description{
Extracts the conditional Stage-2 inverse-Hessian covariance matrix from a
fitted longitudinal functional joint model.
}
\usage{
\method{vcov}{lfjm}(object, ...)
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
Let \eqn{\boldsymbol{\vartheta}} denote the Stage-2 optimization parameter
vector and let
\eqn{Q(\boldsymbol{\vartheta})=-\widetilde{h}(\boldsymbol{\vartheta})} be
the minimized negative adjusted profile h-likelihood. When standard errors
are requested, the fitting procedure evaluates the numerical Hessian at the
optimizer:
\deqn{
\widehat{\mathcal{I}}
=
\left.
\nabla_{\boldsymbol{\vartheta}}^2
Q(\boldsymbol{\vartheta})
\right|_{\boldsymbol{\vartheta}
=\widehat{\boldsymbol{\vartheta}}}.
}{
I_hat = Hessian of Q(vartheta) evaluated at vartheta_hat.
}
The returned covariance matrix is
\deqn{
\widehat{\mathcal{V}}
=
\widehat{\mathcal{I}}^{-1}.
}{
V_hat = inverse(I_hat).
}
The ordinary matrix inverse is used when possible. If it fails, the fitting
procedure attempts to use a generalized inverse. Consequently, an
ill-conditioned fit can yield a covariance matrix that should be interpreted
cautiously even when extraction succeeds.

The matrix is returned on the \emph{optimization scale}. Its row and column
names follow the parameter order
\deqn{
\boldsymbol{\vartheta}
=
\left(
\boldsymbol{\beta}^{\mathsf{T}},
\log\sigma_e,
\log\sigma_u,
\gamma_0,
\gamma_Q,
\alpha,
\log\tau
\right)^{\mathsf{T}},
}{
vartheta = (beta', log_sigma_e, log_sigma_u, gamma0, gamma_Q,
alpha, log_Tau)'.
}
Here, \eqn{\boldsymbol{\beta}} contains, in order, the longitudinal
intercept, visit-time effect, scalar-covariate effect, retained between-level
score coefficients, and retained within-level score coefficients. Their
matrix names are \code{"(Intercept)"}, \code{"obstime"}, \code{"Q1"},
\code{"Bscore1"}, and so forth, followed by \code{"Uscore1"} and so forth.
The remaining names are \code{"log_sigma_e"}, \code{"log_sigma_u"},
\code{"gamma0"}, \code{"gamma_Q"}, \code{"alpha"}, and
\code{"log_Tau"}.

In particular, the entries associated with the residual standard deviation,
random-intercept standard deviation, and Weibull shape parameter do not
directly give natural-scale variances. Let
\deqn{
g(\boldsymbol{\vartheta})
=
\left(
\boldsymbol{\beta}^{\mathsf{T}},
\sigma_e,
\sigma_u,
\gamma_0,
\gamma_Q,
\alpha,
\tau
\right)^{\mathsf{T}}
}{
g(vartheta) = (beta', sigma_e, sigma_u, gamma0, gamma_Q,
alpha, Tau)'
}
be the natural-scale transformation, with
\eqn{\sigma_e=\exp(\log\sigma_e)},
\eqn{\sigma_u=\exp(\log\sigma_u)}, and
\eqn{\tau=\exp(\log\tau)}. If \eqn{\boldsymbol{J}} is the Jacobian of
\eqn{g}, the first-order delta-method covariance matrix is
\deqn{
\widehat{\mathcal{V}}_{\mathrm{natural}}
=
\boldsymbol{J}
\widehat{\mathcal{V}}
\boldsymbol{J}^{\mathsf{T}}.
}{
V_natural = J V_hat J'.
}
The relevant diagonal elements of \eqn{\boldsymbol{J}} are
\eqn{\widehat{\sigma}_e}, \eqn{\widehat{\sigma}_u}, and
\eqn{\widehat{\tau}}; all other diagonal elements are one.
\code{\link{summary.lfjm}} applies these univariate delta-method
transformations when reporting the standard errors of the three positive
parameters.

This covariance matrix is based on the curvature of the adjusted profile
h-likelihood. It is not a sandwich covariance estimator and does not
propagate uncertainty from the Stage-1 estimation of the mean surface,
eigensystems, or functional scores. The matrix therefore describes
conditional Stage-2 uncertainty under the fitted APHL criterion.
}
\value{
A named square numeric matrix with one row and column for every Stage-2
optimization parameter. The matrix is expressed on the optimization scale.
}
\examples{
\dontrun{
simulated_data <- simulate_lfjm(
  setting = "longitudinal",
  number_subjects = 100,
  number_visits = 6,
  seed = 2026
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
  grid = simulated_data$grid,
  control = lfjm_control(
    standard_errors = TRUE
  )
)

covariance_optimizer <- vcov(fit)
sqrt(diag(covariance_optimizer))

# Transform the complete covariance matrix to the natural parameter scale.
jacobian <- diag(nrow(covariance_optimizer))
dimnames(jacobian) <- dimnames(covariance_optimizer)

jacobian["log_sigma_e", "log_sigma_e"] <- fit$sigma_e
jacobian["log_sigma_u", "log_sigma_u"] <- fit$sigma_u
jacobian["log_Tau", "log_Tau"] <- fit$Tau

covariance_natural <-
  jacobian %*%
  covariance_optimizer %*%
  t(jacobian)

natural_names <- rownames(covariance_optimizer)
natural_names[natural_names == "log_sigma_e"] <- "sigma_e"
natural_names[natural_names == "log_sigma_u"] <- "sigma_u"
natural_names[natural_names == "log_Tau"] <- "Tau"

dimnames(covariance_natural) <- list(
  natural_names,
  natural_names
)

sqrt(diag(covariance_natural))
}
}
