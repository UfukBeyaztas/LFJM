\name{print.summary.lfjm}
\alias{print.summary.lfjm}
\title{Print a Summary of an LFJM Fit}
\description{
Prints the model dimensions, adjusted profile h-likelihood value, optimizer
convergence code, and conditional Stage-2 coefficient table stored in an
\code{"summary.lfjm"} object.
}
\usage{
\method{print}{summary.lfjm}(x, digits = 4, ...)
}
\arguments{
  \item{x}{
  An object of class \code{"summary.lfjm"}, typically returned by applying
  \code{\link{summary.lfjm}} to a fitted \code{"lfjm"} object.
  }

  \item{digits}{
  A numeric value specifying the number of significant digits used to display
  the adjusted log-likelihood and coefficient table.
  }

  \item{\dots}{
  Additional arguments. They are currently not used.
  }
}
\details{
The method prints a summary of the two-stage LFPCA--APHL fit. The first part
reports the numbers of subjects and longitudinal records, the retained
between- and within-subject LFPCA dimensions, the maximized adjusted profile
h-likelihood, and the optimizer convergence code.

At Stage 1, longitudinal functional principal component analysis decomposes
the repeatedly observed functional predictor into persistent between-subject
and visit-specific within-subject components. At Stage 2, the retained scores
enter the longitudinal functional joint model as fixed covariates, and the
model parameters are estimated by adjusted profile h-likelihood. For a
random-intercept model, the criterion has the form
\deqn{
\widetilde{h}(\boldsymbol{\psi})
=
\sum_{i=1}^{N}
\left[
h_i\{\widehat{u}_i(\boldsymbol{\psi});\boldsymbol{\psi}\}
- \frac{1}{2}
\log\left\{\frac{H_i(\boldsymbol{\psi})}{2\pi}\right\}
\right],
}{
h_tilde(psi) = sum_i [h_i(u_i_hat(psi); psi)
- 0.5 log\{H_i(psi)/(2 pi)\}].
}
Here, \eqn{\widehat{u}_i(\boldsymbol{\psi})} is the profiled
subject-specific random-intercept mode and \eqn{H_i(\boldsymbol{\psi})} is
the corresponding positive curvature. The displayed value is therefore the
fitted APHL criterion and is not, in general, an exactly marginalized
log-likelihood. The Laplace adjustment is exact when the current-value
association parameter is zero.

The coefficient table is constructed by \code{\link{summary.lfjm}} and
contains the following columns:
\describe{
  \item{\code{Estimate}}{
  Stage-2 parameter estimates. Regression, score, and survival coefficients
  are displayed on their fitted scales. The residual standard deviation
  \code{sigma_e}, random-intercept standard deviation \code{sigma_u}, and
  Weibull shape parameter \code{Tau} are displayed on their natural positive
  scales.
  }
  \item{\code{Std.Error}}{
  Conditional standard errors derived from the inverse numerical Hessian of
  the Stage-2 objective. For \code{sigma_e}, \code{sigma_u}, and \code{Tau},
  the Hessian is formed on the log scale and the reported standard errors are
  transformed to the natural scale by the delta method.
  }
  \item{\code{z.value}}{
  The estimate divided by its conditional standard error.
  }
  \item{\code{p.value}}{
  The two-sided normal-approximation p-value,
  \eqn{2\Phi(-|z|)}.
  }
}

These uncertainty summaries condition on the Stage-1 estimated mean surface,
eigenfunctions, and predicted scores. They do not propagate uncertainty from
the LFPCA stage. If the covariance matrix was not computed, or could not be
obtained, the standard errors, z statistics, and p-values are
\code{NA}. The \code{digits} argument changes only the printed
representation; it does not modify the stored summary object.
}
\value{
The original summary object \code{x}, returned invisibly. The method is
called for its printed side effect.
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

fit_summary <- summary(fit)
print(fit_summary)
print(fit_summary, digits = 3)

invisible_summary <- print(fit_summary)
identical(invisible_summary, fit_summary)
}
}
