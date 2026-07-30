\name{lfjm_control}
\alias{lfjm_control}
\title{Set Control Parameters for LFJM Estimation}
\description{
Constructs a control object containing tuning parameters for the longitudinal
functional principal component analysis and numerical settings for fitting a
longitudinal functional joint model.
}
\usage{
lfjm_control(
  pve_between = 0.95,
  pve_within = 0.95,
  components_between = NULL,
  components_within = NULL,
  mean_degree = 3,
  presmooth = TRUE,
  smoothing_df = 12,
  use_functional_slope = FALSE,
  quadrature_points = 32,
  max_iterations = 500,
  multistart = TRUE,
  standard_errors = TRUE
)
}
\arguments{
  \item{pve_between}{
  A numeric value in \eqn{(0,1]} giving the target proportion of variance
  explained for selecting between-subject functional principal components.
  This setting is used when \code{components_between = NULL}.
  }

  \item{pve_within}{
  A numeric value in \eqn{(0,1]} giving the target proportion of variance
  explained for selecting within-subject functional principal components.
  This setting is used when \code{components_within = NULL}.
  }

  \item{components_between}{
  \code{NULL}, or a positive integer specifying the number of between-subject
  functional principal components to retain. A supplied value overrides
  selection based on \code{pve_between}.
  }

  \item{components_within}{
  \code{NULL}, or a positive integer specifying the number of within-subject
  functional principal components to retain. A supplied value overrides
  selection based on \code{pve_within}.
  }

  \item{mean_degree}{
  A non-negative number controlling the polynomial degree used to estimate
  the functional mean surface over visit time. The effective degree is also
  restricted by the number of distinct visit times and has a maximum of 3.
  }

  \item{presmooth}{
  A logical value indicating whether each observed functional curve should be
  smoothed before the longitudinal functional principal component analysis.
  }

  \item{smoothing_df}{
  The degrees of freedom supplied to \code{\link[stats]{smooth.spline}} when
  \code{presmooth = TRUE}.
  }

  \item{use_functional_slope}{
  A logical value indicating whether the between-subject decomposition should
  use subject-specific functional intercepts and functional slopes. If
  \code{FALSE}, only subject-specific functional intercepts are used.
  }

  \item{quadrature_points}{
  A numeric value of at least 2 giving the number of Gauss--Legendre
  quadrature points used to approximate the cumulative hazard contribution.
  The value is stored as an integer.
  }

  \item{max_iterations}{
  A positive numeric value giving the maximum number of iterations allowed
  for each bounded optimization run. The value is stored as an integer.
  }

  \item{multistart}{
  A logical value indicating whether optimization should be performed from
  multiple initial values for the association parameter. The solution with
  the best objective value is retained and further refined.
  }

  \item{standard_errors}{
  A logical value indicating whether conditional standard errors should be
  calculated from a numerical Hessian after model fitting.
  }
}
\details{
The longitudinal functional joint model is estimated in two stages. In the
first stage, longitudinal functional principal component analysis separates
the functional variation into between-subject and within-subject components.
If a component count is not fixed explicitly, the smallest number of
components whose cumulative eigenvalue proportion reaches the corresponding
\code{pve_between} or \code{pve_within} threshold is retained.

With \code{use_functional_slope = TRUE}, the between-subject process is
represented using a subject-specific functional intercept and slope,
\deqn{
B_i(s;t) = B_i^{(0)}(s) + (t-\bar{t})B_i^{(1)}(s),
}
where \eqn{s} indexes the functional domain and \eqn{t} is visit time. With
the default setting, the between-subject decomposition is based on a
subject-specific functional intercept only.

In the second stage, the joint model is fitted using an adjusted profile
h-likelihood. Gauss--Legendre quadrature approximates the time integral in
the survival contribution. Increasing \code{quadrature_points} may improve
the numerical approximation at the cost of additional computation.
\code{multistart} and \code{max_iterations} affect numerical optimization,
not the statistical specification of the model.

When requested, standard errors are based on the numerical Hessian from the
second-stage fit. They are conditional on the estimated functional quantities
from the first stage and therefore do not propagate first-stage estimation
uncertainty.

Logical controls are stored using \code{isTRUE}; consequently, only an exact
\code{TRUE} value activates the corresponding option.
}
\value{
An object of class \code{"lfjm_control"}, represented by a list with the
following components:
\describe{
  \item{\code{pve_between}}{
  Target proportion of between-subject variance explained.
  }
  \item{\code{pve_within}}{
  Target proportion of within-subject variance explained.
  }
  \item{\code{components_between}}{
  Fixed number of between-subject components, or \code{NULL}.
  }
  \item{\code{components_within}}{
  Fixed number of within-subject components, or \code{NULL}.
  }
  \item{\code{mean_degree}}{
  Control for the polynomial degree of the functional mean surface.
  }
  \item{\code{presmooth}}{
  Whether functional observations are presmoothed.
  }
  \item{\code{smoothing_df}}{
  Degrees of freedom used for presmoothing.
  }
  \item{\code{use_functional_slope}}{
  Whether functional slopes are included in the between-subject
  decomposition.
  }
  \item{\code{quadrature_points}}{
  Integer number of quadrature points.
  }
  \item{\code{max_iterations}}{
  Integer maximum number of optimization iterations.
  }
  \item{\code{multistart}}{
  Whether multiple optimization starts are used.
  }
  \item{\code{standard_errors}}{
  Whether conditional standard errors are requested.
  }
}
}
\examples{
\dontrun{
default_control <- lfjm_control()
default_control

custom_control <- lfjm_control(
  pve_between = 0.99,
  pve_within = 0.95,
  components_between = 3,
  components_within = 2,
  mean_degree = 3,
  presmooth = TRUE,
  smoothing_df = 12,
  use_functional_slope = FALSE,
  quadrature_points = 48,
  max_iterations = 1000,
  multistart = TRUE,
  standard_errors = TRUE
)

custom_control$components_between
custom_control$quadrature_points
}
}
