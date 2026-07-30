\name{lfjm_lfpca}
\alias{lfjm_lfpca}
\title{Longitudinal Functional Principal Component Analysis}
\description{
Performs the first-stage longitudinal functional principal component analysis
used by the longitudinal functional joint model. The procedure estimates a
time-varying functional mean, separates between-subject and within-subject
variation, and obtains weighted-grid eigenfunctions, scores, and reconstructed
functional components.
}
\usage{
lfjm_lfpca(
  functional_data,
  subject,
  time,
  grid,
  pve_between = 0.95,
  pve_within = 0.95,
  components_between = NULL,
  components_within = NULL,
  mean_degree = 4,
  presmooth = TRUE,
  smoothing_df = 12,
  use_functional_slope = FALSE
)
}
\arguments{
  \item{functional_data}{
  A numeric matrix containing the repeatedly observed functional predictor.
  Rows correspond to subject visits and columns correspond to the locations
  in \code{grid}.
  }

  \item{subject}{
  A vector of subject identifiers with one element for each row of
  \code{functional_data}. Repeated identifiers indicate multiple visits from
  the same subject.
  }

  \item{time}{
  A numeric vector of visit times with one element for each row of
  \code{functional_data}.
  }

  \item{grid}{
  A finite, strictly increasing numeric vector giving the common functional
  grid. Its length must equal the number of columns of
  \code{functional_data}.
  }

  \item{pve_between}{
  A numeric value in \eqn{(0,1]} giving the target proportion of
  between-subject variance explained. It is used when
  \code{components_between = NULL}.
  }

  \item{pve_within}{
  A numeric value in \eqn{(0,1]} giving the target proportion of
  within-subject variance explained. It is used when
  \code{components_within = NULL}.
  }

  \item{components_between}{
  \code{NULL}, or a positive integer specifying the number of
  between-subject components to retain. A supplied value overrides
  \code{pve_between}.
  }

  \item{components_within}{
  \code{NULL}, or a positive integer specifying the number of within-subject
  components to retain. A supplied value overrides \code{pve_within}.
  }

  \item{mean_degree}{
  A non-negative number controlling the polynomial degree used for the
  time-varying functional mean. The effective degree cannot exceed 3 or the
  number of distinct visit times minus one.
  }

  \item{presmooth}{
  A logical value indicating whether each observed curve is smoothed before
  the decomposition.
  }

  \item{smoothing_df}{
  The degrees of freedom supplied to \code{\link[stats]{smooth.spline}} when
  \code{presmooth = TRUE}.
  }

  \item{use_functional_slope}{
  A logical value indicating whether the between-subject process should
  contain subject-specific functional intercepts and functional slopes. If
  \code{FALSE}, the between-subject process is represented by a functional
  intercept only.
  }
}
\details{
The working two-level representation of the repeatedly observed functional
predictor is
\deqn{
X_{ij}(s) = \eta(s,t_{ij}) + B_i(s) + U_{ij}(s),
}
where \eqn{\eta(s,t)} is the population mean surface, \eqn{B_i(s)} is the
stable subject-level component, and \eqn{U_{ij}(s)} is the visit-specific
component. This separation permits the later joint model to assign distinct
functional effects to persistent between-subject variation and changing
within-subject variation.

The function first standardizes visit time as
\deqn{
z_{ij} = \{t_{ij}-\bar{t}\}/s_t.
}
If the sample standard deviation of visit time is zero or numerically
negligible, it is replaced by one. At every grid location, the mean surface is
estimated by polynomial regression on standardized time. Its effective degree
is the minimum of 3, \code{mean_degree}, and the number of distinct visit
times minus one. The fitted mean surface is then removed from every curve.

When \code{use_functional_slope = FALSE}, the subject-level curve is estimated
by the average centered curve for that subject. The within-subject curves are
the visit-level residuals obtained after removing this subject average. Thus,
on the observed grid,
\deqn{
\widetilde{X}_{ij}(s) =
\widehat{B}_i^{(0)}(s) + \widehat{U}_{ij}(s).
}

When \code{use_functional_slope = TRUE}, a functional intercept and slope are
estimated separately for each subject by regressing its centered curves on
\eqn{(1,z_{ij})} at all grid locations. The decomposition becomes
\deqn{
\widetilde{X}_{ij}(s) =
\widehat{B}_i^{(0)}(s) +
z_{ij}\widehat{B}_i^{(1)}(s) +
\widehat{U}_{ij}(s).
}
For a subject with only one visit, the observed centered curve is used as the
functional intercept and the functional slope is set to zero.

Separate weighted principal component analyses are applied to the
between-subject and within-subject quantities. If \eqn{\mathbf{W}} is the
diagonal matrix of trapezoidal integration weights, the eigenfunctions are
normalized in the discrete functional inner product
\deqn{
\langle f,g\rangle_{\mathbf{W}} =
f^\mathsf{T}\mathbf{W}g.
}
Consequently, the returned eigenfunctions are orthonormal with respect to the
quadrature weights rather than the unweighted Euclidean inner product. With
functional slopes, the between-subject analysis is performed on the
concatenated intercept and slope curves in the direct-sum space, using the
trapezoidal weights for both parts. Each between eigencomponent therefore
consists of the paired functions returned in \code{phi0} and \code{phi1}.

Unless the number of components is supplied explicitly, the smallest number
whose cumulative eigenvalue proportion reaches the specified proportion of
variance explained is retained. Eigenvalues smaller than
\eqn{10^{-10}} times the largest eigenvalue are treated as numerically zero.
The returned visit-level quantities remain in the original input row order,
whereas subject-level quantities follow the sorted unique subject identifiers
in \code{ids}.
}
\value{
An object of class \code{"lfjm_lfpca"}, represented by a list with the
following components:
\describe{
  \item{\code{call}}{
  The matched function call.
  }
  \item{\code{ids}}{
  The sorted unique subject identifiers.
  }
  \item{\code{subject_index}}{
  Integer indices mapping the input rows to the entries of \code{ids}.
  }
  \item{\code{grid}}{
  The functional grid.
  }
  \item{\code{sgrid}}{
  A duplicate of \code{grid} retained for subsequent model-fitting routines.
  }
  \item{\code{weights}}{
  Trapezoidal integration weights on the functional grid.
  }
  \item{\code{phi0}}{
  A grid-by-\code{KB} matrix containing the functional-intercept part of the
  between-subject eigenfunctions.
  }
  \item{\code{phi1}}{
  A grid-by-\code{KB} matrix containing the functional-slope part of the
  between-subject eigenfunctions. It is zero when
  \code{use_functional_slope = FALSE}.
  }
  \item{\code{phiU}}{
  A grid-by-\code{KU} matrix of within-subject eigenfunctions.
  }
  \item{\code{scoreB_subject}}{
  A subject-by-\code{KB} matrix of between-subject scores. Its row names are
  the subject identifiers.
  }
  \item{\code{scoreU_visit}}{
  A visit-by-\code{KU} matrix of within-subject scores in input row order.
  }
  \item{\code{B0_hat_subject}}{
  A subject-by-grid matrix of reconstructed subject-specific functional
  intercepts.
  }
  \item{\code{B_hat_visit}}{
  A visit-by-grid matrix of reconstructed between-subject curves evaluated at
  each visit time.
  }
  \item{\code{U_hat_visit}}{
  A visit-by-grid matrix of reconstructed within-subject curves.
  }
  \item{\code{eigB}}{
  The retained between-subject eigenvalues.
  }
  \item{\code{eigU}}{
  The retained within-subject eigenvalues.
  }
  \item{\code{KB}}{
  The number of retained between-subject components.
  }
  \item{\code{KU}}{
  The number of retained within-subject components.
  }
  \item{\code{time_mean}}{
  The sample mean of the visit times used for standardization.
  }
  \item{\code{time_sd}}{
  The sample standard deviation used for time standardization, or one when
  the original standard deviation is unavailable or numerically negligible.
  }
  \item{\code{eta_coef}}{
  The polynomial coefficient matrix for the estimated functional mean
  surface.
  }
  \item{\code{eta_degree}}{
  The effective polynomial degree used for the mean surface.
  }
  \item{\code{center_B0}}{
  The gridwise centering curve used in the between functional-intercept
  analysis.
  }
  \item{\code{center_B1}}{
  The gridwise centering curve used in the between functional-slope analysis.
  It is zero when functional slopes are not used.
  }
  \item{\code{center_U}}{
  The gridwise centering curve used in the within-subject analysis.
  }
  \item{\code{presmooth}}{
  The supplied curve-presmoothing indicator.
  }
  \item{\code{smooth_df}}{
  The supplied smoothing degrees of freedom.
  }
  \item{\code{use_slope}}{
  The supplied functional-slope indicator.
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

print(lfpca_fit)
c(between = lfpca_fit$KB, within = lfpca_fit$KU)
dim(lfpca_fit$scoreB_subject)
dim(lfpca_fit$scoreU_visit)
}
}
