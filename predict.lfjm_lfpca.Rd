\name{predict.lfjm_lfpca}
\alias{predict.lfjm_lfpca}
\title{Project Functional Histories onto an LFPCA Fit}
\description{
Projects new repeatedly observed functional histories onto the mean surface
and between- and within-level eigensystems of a fitted longitudinal functional
principal component analysis.
}
\usage{
\method{predict}{lfjm_lfpca}(
  object,
  functional_data,
  subject,
  time,
  ...
)
}
\arguments{
  \item{object}{
  An object of class \code{"lfjm_lfpca"}, typically returned by
  \code{\link{lfjm_lfpca}}.
  }

  \item{functional_data}{
  A numeric matrix containing the new functional observations. Rows correspond
  to subject visits and columns must correspond, in the same order, to the
  fitted grid stored in \code{object$grid}.
  }

  \item{subject}{
  A vector of subject identifiers with one element for each row of
  \code{functional_data}. Repeated identifiers define the functional history
  supplied for each subject.
  }

  \item{time}{
  A numeric vector of visit times with one element for each row of
  \code{functional_data}.
  }

  \item{\dots}{
  Additional arguments. They are currently not used.
  }
}
\details{
This method computes functional principal component scores using an existing
LFPCA fit. It does not re-estimate the functional mean surface, centering
curves, eigenfunctions, eigenvalues, or truncation dimensions.

The new functional observations must be evaluated at the same grid locations
and in the same column order as the training curves. The method checks the
number of grid columns but, because a separate new grid is not supplied, it
cannot verify the locations themselves. If the original analysis used curve
presmoothing, each new curve is smoothed using the fitted grid and the stored
degrees of freedom before projection.

Visit times are standardized using the training-sample quantities,
\deqn{
z_{ij}^{\mathrm{new}} =
\frac{t_{ij}^{\mathrm{new}}-\overline{t}_{\mathrm{fit}}}
{s_{t,\mathrm{fit}}}.
}{
z_ij_new = (t_ij_new - time_mean_fit) / time_sd_fit.
}
The fitted polynomial mean surface is evaluated at these standardized times
and removed from the new curves:
\deqn{
\widetilde{X}_{ij}^{\mathrm{new}}(s)
=
X_{ij}^{\mathrm{new}}(s)
-\widehat{\eta}(s,t_{ij}^{\mathrm{new}}).
}{
X_tilde_ij_new(s) =
X_ij_new(s) - eta_hat(s, t_ij_new).
}

\subsection{Functional-intercept projection}{
When the fitted object has \code{use_slope = FALSE}, the subject-level
component is the average of the centered curves supplied for that subject,
\deqn{
\widehat{B}^{(0)}_i(s)
=
\frac{1}{n_i}\sum_{j=1}^{n_i}
\widetilde{X}_{ij}^{\mathrm{new}}(s),
}{
B0_i_hat(s) = mean_j X_tilde_ij_new(s),
}
and the visit-level component is
\deqn{
\widehat{U}_{ij}(s)
=
\widetilde{X}_{ij}^{\mathrm{new}}(s)
-\widehat{B}^{(0)}_i(s).
}{
U_ij_hat(s) = X_tilde_ij_new(s) - B0_i_hat(s).
}

Let \eqn{\mathbf W} be the diagonal matrix of trapezoidal weights, let
\eqn{\boldsymbol{\Phi}_0} contain the fitted between-level eigenfunctions,
and let \eqn{\boldsymbol{\Phi}_U} contain the fitted within-level
eigenfunctions. The projected scores are
\deqn{
\widehat{\boldsymbol{\xi}}_i
=
\{\widehat{\mathbf B}^{(0)}_i-\mathbf c_{B0}\}^{\mathsf T}
\mathbf W\boldsymbol{\Phi}_0,
\qquad
\widehat{\boldsymbol{\zeta}}_{ij}
=
\{\widehat{\mathbf U}_{ij}-\mathbf c_U\}^{\mathsf T}
\mathbf W\boldsymbol{\Phi}_U,
}{
xi_i_hat = (B0_i_hat - center_B0)' W Phi0,
zeta_ij_hat = (U_ij_hat - center_U)' W PhiU.
}
The centering curves \eqn{\mathbf c_{B0}} and \eqn{\mathbf c_U} are those
estimated in the original fit.
}

\subsection{Functional-intercept and slope projection}{
When the fitted object has \code{use_slope = TRUE}, the centered curves of
each new subject are regressed pointwise on \eqn{(1,z_{ij}^{\mathrm{new}})}.
This produces a subject-specific functional intercept
\eqn{\widehat{B}^{(0)}_i}, a functional slope
\eqn{\widehat{B}^{(1)}_i}, and within-subject residual curves. For a subject
with one supplied visit, the centered curve is used as the functional
intercept and the functional slope is set to zero.

The between-level score projection is performed in the direct-sum
intercept--slope space:
\deqn{
\widehat{\boldsymbol{\xi}}_i
=
\left[
\begin{array}{c}
\widehat{\mathbf B}^{(0)}_i-\mathbf c_{B0}\\
\widehat{\mathbf B}^{(1)}_i-\mathbf c_{B1}
\end{array}
\right]^{\mathsf T}
\left(\mathbf W\mathbin{\oplus}\mathbf W\right)
\left[
\begin{array}{c}
\boldsymbol{\Phi}_0\\
\boldsymbol{\Phi}_1
\end{array}
\right].
}{
xi_i_hat =
c(B0_i_hat - center_B0, B1_i_hat - center_B1)'
diag(W, W) rbind(Phi0, Phi1).
}
The within-level scores use the same weighted projection as in the
functional-intercept case, applied to the residual curves from the
subject-specific intercept--slope regression.
}

The between-level scores depend on the complete set of visits supplied for a
new subject. Consequently, adding another visit for that subject can change
both its estimated between-level score and the within-level scores of its
earlier visits. Subject-level output follows the sorted unique identifiers,
whereas visit-level output remains in the input row order.
}
\value{
A list with the following components:
\describe{
  \item{\code{ids}}{
  The sorted unique subject identifiers in the supplied data.
  }
  \item{\code{subject_index}}{
  Integer indices mapping each input row to the corresponding entry in
  \code{ids}.
  }
  \item{\code{scoreB_subject}}{
  A subject-by-\code{object$KB} matrix of projected between-level scores.
  Row names are the subject identifiers.
  }
  \item{\code{scoreU_visit}}{
  A visit-by-\code{object$KU} matrix of projected within-level scores in the
  original input row order.
  }
}
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

new_subject <- rep(c("A", "B"), each = 3)
new_time <- rep(c(0, 1, 2), times = 2)
new_subject_effect <- rep(c(-0.30, 0.50), each = 3)

new_functional_data <-
  matrix(
    0.25 * sin(pi * grid),
    nrow = length(new_subject),
    ncol = length(grid),
    byrow = TRUE
  ) +
  outer(new_subject_effect, sin(pi * grid)) +
  outer(0.20 * new_time, cos(2 * pi * grid)) +
  matrix(
    rnorm(length(new_subject) * length(grid), sd = 0.05),
    nrow = length(new_subject)
  )

projected_scores <- predict(
  lfpca_fit,
  functional_data = new_functional_data,
  subject = new_subject,
  time = new_time
)

projected_scores$ids
projected_scores$subject_index
projected_scores$scoreB_subject
projected_scores$scoreU_visit
}
}
