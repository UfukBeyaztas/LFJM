\name{print.lfjm_lfpca}
\alias{print.lfjm_lfpca}
\title{Print a Longitudinal Functional PCA Fit}
\description{
Prints a compact overview of a longitudinal functional principal component
analysis, including the sample size, grid size, retained component counts, and
curve-presmoothing status.
}
\usage{
\method{print}{lfjm_lfpca}(x, ...)
}
\arguments{
  \item{x}{
  An object of class \code{"lfjm_lfpca"}, typically returned by
  \code{\link{lfjm_lfpca}}.
  }

  \item{\dots}{
  Additional arguments. They are currently not used.
  }
}
\details{
Longitudinal functional principal component analysis separates the centered
functional predictor into a persistent subject-level process and a
visit-specific process,
\deqn{
\widetilde{X}_{ij}(s) =
B_i(s)+U_{ij}(s).
}{
X_tilde_ij(s) = B_i(s) + U_ij(s).
}
The two processes are represented using separate truncated eigensystems. The
printed between- and within-level component counts are the retained dimensions
\eqn{K_B} and \eqn{K_U}, stored in \code{x$KB} and \code{x$KU}.

The method prints the following fields:
\describe{
  \item{\code{Subjects}}{
  The number of unique subjects in the fitted analysis, calculated as
  \code{length(x$ids)}.
  }
  \item{\code{Functional grid size}}{
  The number of locations in the common functional grid, calculated as
  \code{length(x$grid)}.
  }
  \item{\code{Between-level components}}{
  The number of retained subject-level functional principal components,
  \code{x$KB}.
  }
  \item{\code{Within-level components}}{
  The number of retained visit-level functional principal components,
  \code{x$KU}.
  }
  \item{\code{Curve presmoothing}}{
  \code{"yes"} when \code{x$presmooth} is true and \code{"no"} otherwise.
  Presmoothing is applied to each observed functional curve before estimation
  of the mean surface and the two-level decomposition.
  }
}

When the LFPCA was fitted with functional slopes, each between-level
eigencomponent is a paired functional intercept and functional slope in a
direct-sum space. In that case, the displayed value \code{x$KB} counts these
paired components, rather than counting the intercept and slope functions
separately.
}
\value{
The original fitted object \code{x}, returned invisibly. The method is called
for its printed side effect.
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

invisible_fit <- print(lfpca_fit)
identical(invisible_fit, lfpca_fit)
}
}
