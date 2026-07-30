\name{print.lfjm_control}
\alias{print.lfjm_control}
\title{Print LFJM Control Parameters}
\description{
Prints a compact selection of the component-selection, mean-surface,
quadrature, and optimization settings stored in an LFJM control object.
}
\usage{
\method{print}{lfjm_control}(x, ...)
}
\arguments{
  \item{x}{
  An object of class \code{"lfjm_control"}, created by
  \code{\link{lfjm_control}}.
  }

  \item{\dots}{
  Additional arguments. They are currently not used.
  }
}
\details{
The printed output summarizes settings governing both stages of LFJM
estimation. The between- and within-level proportion-of-variance-explained
thresholds are printed from \code{x$pve_between} and
\code{x$pve_within}. When the corresponding component count is \code{NULL},
the method prints \code{"selected by PVE"}. In that case, Stage 1 retains the
smallest number of components satisfying
\deqn{
\frac{\sum_{k=1}^{K}\widehat{\lambda}_k}
{\sum_{k}\widehat{\lambda}_k}
\geq p,
}{
sum of the first K eigenvalues / sum of all eigenvalues >= p,
}
where \eqn{p} is the printed PVE threshold. If a component count is supplied
explicitly, that value is printed and takes precedence over the corresponding
PVE rule.

\code{Mean degree} is the control value stored in \code{x$mean_degree} for the
polynomial time trend in the functional mean surface. The printed value is the
requested control, not necessarily the effective fitted degree. The
mean-surface estimator also limits the degree to 3 and to the number of
distinct visit times minus one.

\code{Quadrature points} gives the number of Gauss--Legendre nodes used to
approximate the cumulative-hazard integral in the Stage-2 adjusted profile
h-likelihood. Increasing this value can improve the numerical approximation
but increases computational cost. \code{Maximum iterations} is the iteration
limit supplied to each bounded optimization run.

The method prints the following fields:
\describe{
  \item{\code{Between-level PVE}}{
  The target proportion of between-subject functional variance explained.
  }
  \item{\code{Within-level PVE}}{
  The target proportion of within-subject functional variance explained.
  }
  \item{\code{Between components}}{
  The fixed between-level component count, or an indication that PVE-based
  selection is used.
  }
  \item{\code{Within components}}{
  The fixed within-level component count, or an indication that PVE-based
  selection is used.
  }
  \item{\code{Mean degree}}{
  The requested polynomial-degree control for the functional mean surface.
  }
  \item{\code{Quadrature points}}{
  The number of numerical-integration points used in Stage 2.
  }
  \item{\code{Maximum iterations}}{
  The maximum number of iterations for each bounded optimization run.
  }
}
}
\value{
The original control object \code{x}, returned invisibly. The method is called
for its printed side effect.
}
\examples{
\dontrun{
default_control <- lfjm_control()
print(default_control)

fixed_components <- lfjm_control(
  pve_between = 0.99,
  pve_within = 0.95,
  components_between = 3,
  components_within = 2,
  mean_degree = 3,
  quadrature_points = 48,
  max_iterations = 1000
)

print(fixed_components)

invisible_control <- print(fixed_components)
identical(invisible_control, fixed_components)

fixed_components$presmooth
fixed_components$standard_errors
}
}
