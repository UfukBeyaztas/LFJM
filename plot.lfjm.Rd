\name{plot.lfjm}
\alias{plot.lfjm}
\title{Plot Functional Components from an LFJM Fit}
\description{
Produces plots of the estimated functional coefficient functions or the
between- and within-level eigensystems from a fitted longitudinal functional
joint model.
}
\usage{
\method{plot}{lfjm}(
  x,
  type = c("coefficients", "eigenfunctions", "eigenvalues"),
  level = c("between", "within"),
  ...
)
}
\arguments{
  \item{x}{
  An object of class \code{"lfjm"}, typically returned by
  \code{\link{lfjm}}.
  }

  \item{type}{
  A character string specifying the plot to produce. The available choices
  are \code{"coefficients"}, \code{"eigenfunctions"}, and
  \code{"eigenvalues"}. Partial matching is allowed.
  }

  \item{level}{
  A character string specifying whether the between-subject or within-subject
  eigensystem is plotted. The available choices are \code{"between"} and
  \code{"within"}. This argument is validated but otherwise ignored when
  \code{type = "coefficients"}.
  }

  \item{\dots}{
  Additional graphical arguments passed to \code{\link[graphics]{matplot}} for
  coefficient and eigenfunction plots, or to \code{\link[graphics]{plot}} for
  eigenvalue plots. Arguments that duplicate graphical parameters fixed
  internally by the method cannot be used to override those parameters.
  }
}
\details{
The fitted model separates the centered longitudinal functional predictor into
a persistent between-subject process \eqn{B_i(s)} and a visit-specific
within-subject process \eqn{U_{ij}(s)}. Their effects in the longitudinal
submodel are represented by two coefficient functions,
\eqn{\beta_B(s)} and \eqn{\beta_U(s)}. The shared current-value association
routes the between-subject functional effect to the survival submodel through
\eqn{\alpha\beta_B(s)}.

With \code{type = "coefficients"}, the method plots the three estimated
functions stored in the fitted object:
\deqn{
\widehat{\beta}_B(s),\qquad
\widehat{\beta}_U(s),\qquad
\widehat{\alpha}\widehat{\beta}_B(s).
}{
beta_B_hat(s), beta_U_hat(s), and alpha_hat beta_B_hat(s).
}
They are evaluated at \code{x$grid}. The between-subject longitudinal
coefficient is shown by a solid black line, the within-subject longitudinal
coefficient by a dashed blue line, and the induced survival coefficient by a
dotted red line.

The functional coefficients are obtained by mapping the estimated Stage-2
score coefficients back to the functional grid. For the retained centered
representations
\eqn{\widehat{B}_{i,K_B}=\sum_k\widehat{\xi}_{ik}\widehat{\phi}^B_k} and
\eqn{\widehat{U}_{ij,K_U}=\sum_l\widehat{\zeta}_{ijl}\widehat{\phi}^U_l},
the trapezoidal-grid inner product gives
\deqn{
\langle \widehat{B}_{i,K_B},\widehat{\beta}_B\rangle_{\mathbf W}
=
\widehat{\boldsymbol{\xi}}_i^\mathsf{T}
\widehat{\boldsymbol{\theta}}_B,
\qquad
\langle \widehat{U}_{ij,K_U},\widehat{\beta}_U\rangle_{\mathbf W}
=
\widehat{\boldsymbol{\zeta}}_{ij}^\mathsf{T}
\widehat{\boldsymbol{\theta}}_U.
}{
<B_iK_hat, beta_B_hat>_W = xi_i_hat' theta_B_hat,
<U_ijK_hat, beta_U_hat>_W = zeta_ij_hat' theta_U_hat.
}
Thus, the coefficient plot displays the estimated functions themselves; a
subject- or visit-specific functional contribution additionally depends on
the corresponding functional curve or scores.

With \code{type = "eigenfunctions"} or \code{type = "eigenvalues"}, plotting
is delegated to the Stage-1 object stored in \code{x$lfpca}. For
\code{level = "between"}, the eigenfunction plot displays the columns of
\code{x$lfpca$phi0}, namely the functional-intercept part of the retained
between-subject eigenfunctions. If functional slopes were used in Stage 1,
their \code{phi1} parts are not displayed by this method. For
\code{level = "within"}, the columns of \code{x$lfpca$phiU} are displayed.

The eigenvalue plot shows the retained eigenvalues in component order, using
\code{x$lfpca$eigB} or \code{x$lfpca$eigU}. Eigenvalues quantify the variation
associated with the retained directions at the selected level. If no
components were retained at that level, an empty plotting region is produced
with an explanatory title.
}
\value{
The original fitted object \code{x}, returned invisibly. The function is
called for its graphical side effect.
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

plot(fit)
plot(fit, type = "coefficients", ylim = c(-2, 2))
plot(fit, type = "eigenfunctions", level = "between")
plot(fit, type = "eigenfunctions", level = "within")
plot(fit, type = "eigenvalues", level = "between")
plot(fit, type = "eigenvalues", level = "within")
}
}
