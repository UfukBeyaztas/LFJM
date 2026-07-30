\name{plot.lfjm_lfpca}
\alias{plot.lfjm_lfpca}
\title{Plot Longitudinal Functional PCA Eigensystems}
\description{
Plots the retained eigenfunctions or eigenvalues from the between-subject or
within-subject level of a longitudinal functional principal component
analysis.
}
\usage{
\method{plot}{lfjm_lfpca}(
  x,
  type = c("eigenfunctions", "eigenvalues"),
  level = c("between", "within"),
  ...
)
}
\arguments{
  \item{x}{
  An object of class \code{"lfjm_lfpca"}, typically returned by
  \code{\link{lfjm_lfpca}}.
  }

  \item{type}{
  A character string specifying whether to plot the retained
  \code{"eigenfunctions"} or \code{"eigenvalues"}. Partial matching is
  allowed.
  }

  \item{level}{
  A character string specifying the decomposition level. Use
  \code{"between"} for the subject-level eigensystem or \code{"within"} for
  the visit-level eigensystem. Partial matching is allowed.
  }

  \item{\dots}{
  Additional graphical arguments passed to \code{\link[graphics]{matplot}} for
  eigenfunction plots or to \code{\link[graphics]{plot}} for eigenvalue plots.
  Arguments that duplicate graphical parameters fixed internally by the
  method cannot be used to override those parameters.
  }
}
\details{
The longitudinal functional principal component analysis separates the
centered functional predictor into between-subject and within-subject
processes. In the functional-intercept formulation,
\deqn{
\widetilde{X}_{ij}(s) =
B_i(s) + U_{ij}(s),
}{
X_tilde_ij(s) = B_i(s) + U_ij(s),
}
where \eqn{B_i(s)} describes persistent variation among subjects and
\eqn{U_{ij}(s)} describes variation of visit \eqn{j} around the corresponding
subject-level curve.

The two processes are represented by separate truncated Karhunen--Loeve
expansions,
\deqn{
B_i(s) \simeq
\sum_{k=1}^{K_B}\xi_{ik}\phi_k^B(s),
\qquad
U_{ij}(s) \simeq
\sum_{l=1}^{K_U}\zeta_{ijl}\phi_l^U(s).
}{
B_i(s) = sum_k xi_ik phi_k_B(s),
U_ij(s) = sum_l zeta_ijl phi_l_U(s).
}
The between-subject eigenfunctions are stored in \code{x$phi0}, the
within-subject eigenfunctions in \code{x$phiU}, and their retained eigenvalues
in \code{x$eigB} and \code{x$eigU}, respectively.

With \code{type = "eigenfunctions"}, all retained eigenfunctions at the
selected level are drawn as lines over \code{x$grid}. The eigenfunctions were
obtained from a weighted-grid decomposition using the trapezoidal integration
weights \eqn{w_d}. Their discrete orthonormality is therefore
\deqn{
\sum_d w_d\phi_k(s_d)\phi_{k'}(s_d)
=
\delta_{kk'},
}{
sum_d w_d phi_k(s_d) phi_k'(s_d) = delta_kk'.
}
They need not be orthonormal under the ordinary unweighted Euclidean inner
product.

For \code{level = "between"}, the method plots \code{x$phi0}. When the LFPCA
was fitted with \code{use_functional_slope = TRUE}, each between-subject
component consists of a functional-intercept part in \code{x$phi0} and a
functional-slope part in \code{x$phi1}; this plotting method displays only the
functional-intercept part. For \code{level = "within"}, it plots
\code{x$phiU}. In the functional-slope case, the between-level eigenvalues
refer to the complete intercept--slope components in the direct-sum space,
even though only their intercept parts are drawn in the eigenfunction plot.

With \code{type = "eigenvalues"}, the method plots the retained eigenvalues
against their component indices. Eigenvalues quantify the level-specific
variation associated with each retained direction and appear in decreasing
order. The plot contains only the components retained by the fixed-component
or proportion-of-variance rule used in \code{\link{lfjm_lfpca}}; it does not
display discarded eigenvalues.

If no components are available at the requested level, the method creates an
empty plotting region with a title stating that no components were retained.
}
\value{
The original \code{"lfjm_lfpca"} object \code{x}, returned invisibly. The
function is called for its graphical side effect.
}
\examples{
\dontrun{
grid <- seq(0, 1, length.out = 31)
subject <- rep(1:6, each = 3)
time <- rep(c(0, 1, 2), times = 6)

subject_effect <- seq(-0.75, 0.75, length.out = 6)[subject]

functional_data <-
  matrix(
    0.25 * sin(pi * grid),
    nrow = length(subject),
    ncol = length(grid),
    byrow = TRUE
  ) +
  outer(subject_effect, sin(pi * grid)) +
  outer(0.20 * time, cos(2 * pi * grid)) +
  matrix(
    rnorm(length(subject) * length(grid), sd = 0.05),
    nrow = length(subject)
  )

lfpca_fit <- lfjm_lfpca(
  functional_data = functional_data,
  subject = subject,
  time = time,
  grid = grid,
  components_between = 2,
  components_within = 2,
  mean_degree = 3,
  presmooth = FALSE
)

plot(lfpca_fit, type = "eigenfunctions", level = "between")
plot(lfpca_fit, type = "eigenfunctions", level = "within")
plot(lfpca_fit, type = "eigenvalues", level = "between")
plot(lfpca_fit, type = "eigenvalues", level = "within")
}
}
