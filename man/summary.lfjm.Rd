\name{summary.lfjm}
\alias{summary.lfjm}
\title{Summarize a Longitudinal Functional Joint Model Fit}
\description{
Constructs a summary of a fitted longitudinal functional joint model,
including conditional Stage-2 coefficient estimates, standard errors, Wald
statistics, fit information, and retained functional component counts.
}
\usage{
\method{summary}{lfjm}(object, ...)
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
The longitudinal functional joint model is fitted in two stages. Stage 1
estimates the mean surface, between- and within-subject eigensystems, and
functional principal component scores. Stage 2 treats these estimated
quantities as fixed covariates and estimates the joint-model parameters by
adjusted profile h-likelihood.

For retained between-subject scores
\eqn{\widehat{\boldsymbol{\xi}}_i} and within-subject scores
\eqn{\widehat{\boldsymbol{\zeta}}_{ij}}, the working longitudinal submodel is
\deqn{
Y_{ij}
=
\beta_0+\beta_{\mathrm{time}}t_{ij}+\beta_Q Q_i
+\widehat{\boldsymbol{\xi}}_i^\mathsf{T}\boldsymbol{\theta}_B
+\widehat{\boldsymbol{\zeta}}_{ij}^\mathsf{T}\boldsymbol{\theta}_U
+u_i+\varepsilon_{ij}.
}{
Y_ij = beta_0 + beta_time t_ij + beta_Q Q_i
+ xi_hat_i' theta_B + zeta_hat_ij' theta_U + u_i + epsilon_ij.
}
The survival hazard is
\deqn{
h_i(t)
=
\tau t^{\tau-1}
\exp\{\gamma_0+\gamma_Q Q_i+\alpha m_i(t)\},
}{
h_i(t) = tau t^(tau-1)
exp\{gamma_0 + gamma_Q Q_i + alpha m_i(t)\},
}
where the shared current value is
\deqn{
m_i(t)
=
\beta_0+\beta_{\mathrm{time}}t+\beta_Q Q_i
+\widehat{\boldsymbol{\xi}}_i^\mathsf{T}\boldsymbol{\theta}_B+u_i.
}{
m_i(t) = beta_0 + beta_time t + beta_Q Q_i
+ xi_hat_i' theta_B + u_i.
}

The rows of the coefficient table are ordered as follows:
\enumerate{
  \item the longitudinal intercept, visit-time coefficient, and scalar
  longitudinal-covariate coefficient, named \code{"(Intercept)"},
  \code{"obstime"}, and \code{"Q1"};
  \item the retained between-level score coefficients, named
  \code{"Bscore1"}, \code{"Bscore2"}, and so forth;
  \item the retained within-level score coefficients, named
  \code{"Uscore1"}, \code{"Uscore2"}, and so forth; and
  \item \code{"sigma_e"}, \code{"sigma_u"}, \code{"gamma0"},
  \code{"gamma_Q"}, \code{"alpha"}, and \code{"Tau"}.
}
Thus, the table reports the finite-dimensional score coefficients
\eqn{\boldsymbol{\theta}_B} and \eqn{\boldsymbol{\theta}_U}. It does not
contain the reconstructed grid-level coefficient functions
\eqn{\widehat{\beta}_B(s)}, \eqn{\widehat{\beta}_U(s)}, or the induced
survival functional effect
\eqn{\widehat{\alpha}\widehat{\beta}_B(s)}. These can be extracted with
\code{\link{coef.lfjm}}.

When \code{object$covariance_matrix} is available, standard errors are
calculated from the diagonal of the inverse numerical Hessian of the
minimized negative Stage-2 criterion. The optimization parameterization uses
\eqn{\log(\sigma_e)}, \eqn{\log(\sigma_u)}, and \eqn{\log(\tau)}. The
corresponding natural-scale standard errors are therefore obtained by the
first-order delta method:
\deqn{
\operatorname{se}(\widehat{\sigma}_e)
=
\widehat{\sigma}_e
\operatorname{se}\{\log(\widehat{\sigma}_e)\},
\qquad
\operatorname{se}(\widehat{\sigma}_u)
=
\widehat{\sigma}_u
\operatorname{se}\{\log(\widehat{\sigma}_u)\},
}{
se(sigma_e_hat) = sigma_e_hat se(log(sigma_e_hat)),
se(sigma_u_hat) = sigma_u_hat se(log(sigma_u_hat)),
}
and
\deqn{
\operatorname{se}(\widehat{\tau})
=
\widehat{\tau}
\operatorname{se}\{\log(\widehat{\tau})\}.
}{
se(tau_hat) = tau_hat se(log(tau_hat)).
}
The remaining standard errors are taken directly from the optimization-scale
covariance matrix. For every displayed parameter estimate
\eqn{\widehat{\vartheta}}, the method forms
\deqn{
z
=
\frac{\widehat{\vartheta}}
{\operatorname{se}(\widehat{\vartheta})},
\qquad
p=2\Phi(-|z|).
}{
z = estimate / standard error,
p = 2 Phi(-abs(z)).
}
If the covariance matrix was not requested or could not be computed, all
standard errors, z statistics, and p-values are returned as \code{NA}.

The reported \code{logLik} is the maximized adjusted profile h-likelihood,
\deqn{
\widetilde{h}(\widehat{\boldsymbol{\psi}})
=
\sum_{i=1}^{N}
\left[
h_i(\widehat{u}_i;\widehat{\boldsymbol{\psi}})
-\frac{1}{2}
\log\left\{\frac{H_i}{2\pi}\right\}
\right],
}{
h_tilde(psi_hat) =
sum_i [h_i(u_i_hat; psi_hat) - 0.5 log\{H_i/(2 pi)\}].
}
Here, \eqn{\widehat{u}_i} is the profiled random-intercept mode and
\eqn{H_i} is its positive curvature. The value is a Laplace-adjusted profile
criterion rather than, in general, an exactly marginalized log-likelihood;
the adjustment is exact when the association parameter is zero.
}
\value{
An object of class \code{"summary.lfjm"}, represented by a list with the
following components:
\describe{
  \item{\code{call}}{
  The matched call used to fit the model.
  }

  \item{\code{coefficients}}{
  A data frame with rows corresponding to the Stage-2 model parameters and
  columns \code{Estimate}, \code{Std.Error}, \code{z.value}, and
  \code{p.value}.
  }

  \item{\code{logLik}}{
  The maximized adjusted profile h-likelihood value.
  }

  \item{\code{objective}}{
  The minimized negative adjusted profile h-likelihood value returned by the
  optimizer.
  }

  \item{\code{convergence}}{
  The optimizer convergence code. Under the optimization routines used by
  \code{\link{lfjm}}, zero ordinarily indicates successful convergence.
  }

  \item{\code{convergence_message}}{
  The optional optimizer convergence message.
  }

  \item{\code{number_subjects}}{
  The number of subject-level survival records.
  }

  \item{\code{number_longitudinal_records}}{
  The number of retained longitudinal observations.
  }

  \item{\code{components_between}, \code{components_within}}{
  The numbers of retained between- and within-subject LFPCA components.
  }
}
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
  grid = simulated_data$grid
)

fit_summary <- summary(fit)
print(fit_summary)

fit_summary$coefficients
fit_summary$convergence
fit_summary$components_between
fit_summary$components_within
}
}
