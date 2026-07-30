\name{print.lfjm}
\alias{print.lfjm}
\title{Print a Longitudinal Functional Joint Model Fit}
\description{
Prints a compact overview of a fitted longitudinal functional joint model,
including sample sizes, retained functional component counts, the adjusted
profile h-likelihood value, and the optimizer convergence code.
}
\usage{
\method{print}{lfjm}(x, ...)
}
\arguments{
  \item{x}{
  An object of class \code{"lfjm"}, typically returned by
  \code{\link{lfjm}}.
  }

  \item{\dots}{
  Additional arguments. They are currently not used.
  }
}
\details{
The method identifies the fitted procedure as LFPCA--APHL, reflecting its
two-stage construction. In Stage 1, longitudinal functional principal
component analysis separates the repeatedly observed functional predictor
into between-subject and within-subject components. The displayed component
counts are
\code{x$lfpca$KB} and \code{x$lfpca$KU}, the numbers of retained directions
at these two levels.

Stage 2 fits the longitudinal and survival submodels through the adjusted
profile h-likelihood. In the random-intercept implementation, the displayed
criterion is
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
\eqn{H_i} is the corresponding positive curvature. This value is stored in
\code{x$loglik} and printed to three decimal places. It is the fitted APHL
criterion, not generally an exactly marginalized log-likelihood; the Laplace
adjustment is exact when the association parameter is zero.

The output contains:
\describe{
  \item{\code{Subjects}}{
  The number of rows in \code{x$survival_data}. The fitted data structure has
  one survival record per subject.
  }
  \item{\code{Longitudinal records}}{
  The number of rows in \code{x$longitudinal_data}.
  }
  \item{\code{Between-level components}}{
  The retained between-subject LFPCA dimension, \code{x$lfpca$KB}.
  }
  \item{\code{Within-level components}}{
  The retained within-subject LFPCA dimension, \code{x$lfpca$KU}.
  }
  \item{\code{APHL}}{
  The maximized adjusted profile h-likelihood stored in \code{x$loglik}.
  }
  \item{\code{Optimizer convergence code}}{
  The code stored in \code{x$convergence}. Under the optimization routines
  used by \code{\link{lfjm}}, a value of zero ordinarily indicates successful
  convergence. A nonzero value should be investigated together with
  \code{x$convergence_message}.
  }
}
}
\value{
The original fitted object \code{x}, returned invisibly. The method is called
for its printed side effect.
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

print(fit)

invisible_fit <- print(fit)
identical(invisible_fit, fit)
}
}
