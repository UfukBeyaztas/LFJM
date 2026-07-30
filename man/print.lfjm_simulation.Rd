\name{print.lfjm_simulation}
\alias{print.lfjm_simulation}
\title{Print Simulated LFJM Data}
\description{
Prints a compact overview of a simulated longitudinal functional joint-model
data set, including its simulation setting, sample dimensions, functional grid
size, and observed event rate.
}
\usage{
\method{print}{lfjm_simulation}(x, ...)
}
\arguments{
  \item{x}{
  An object of class \code{"lfjm_simulation"}, typically returned by
  \code{\link{simulate_lfjm}}.
  }

  \item{\dots}{
  Additional arguments. They are currently not used.
  }
}
\details{
The simulation generator provides two functional-predictor settings. In the
\code{"longitudinal"} setting, the functional predictor contains both a
persistent between-subject process and visit-specific within-subject
variation:
\deqn{
X_{ij}(s) =
\eta(s,t_{ij})+B_i(s)+U_{ij}(s)+e_{ij}(s).
}{
X_ij(s) = eta(s, t_ij) + B_i(s) + U_ij(s) + e_ij(s).
}
In the \code{"baseline"} setting, each subject's first retained functional
curve is repeated over its visits and the within-subject functional component
and its coefficient are set to zero. The \code{Setting} line prints the label
stored in \code{x$case_type}.

The method prints the following fields:
\describe{
  \item{\code{Setting}}{
  The simulation case stored in \code{x$case_type}, currently
  \code{"longitudinal"} or \code{"baseline"}.
  }
  \item{\code{Subjects}}{
  The number of rows in \code{x$data.surv}. The simulation creates one
  survival record per subject.
  }
  \item{\code{Longitudinal records}}{
  The number of rows in \code{x$data.long}. Scheduled visits occurring after a
  subject's observed follow-up time are removed, so this number can be smaller
  than the number of subjects multiplied by the requested number of visits.
  }
  \item{\code{Functional grid size}}{
  The number of locations in \code{x$grid}.
  }
  \item{\code{Event rate}}{
  The arithmetic mean of the binary event indicator in
  \code{x$data.surv$event}, printed to three decimal places:
  \deqn{
  \widehat{p}_{\mathrm{event}}
  =
  \frac{1}{N}\sum_{i=1}^{N}\delta_i.
  }{
  event rate = sum_i event_i / N.
  }
  This is the realized proportion of observed failures rather than the
  uncensored failure-time distribution or a hazard estimate.
  }
}
}
\value{
The original simulation object \code{x}, returned invisibly. The method is
called for its printed side effect.
}
\examples{
\dontrun{
longitudinal_data <- simulate_lfjm(
  setting = "longitudinal",
  number_subjects = 100,
  number_visits = 6
)

print(longitudinal_data)

baseline_data <- simulate_lfjm(
  setting = "baseline",
  number_subjects = 100,
  number_visits = 6
)

print(baseline_data)

invisible_data <- print(longitudinal_data)
identical(invisible_data, longitudinal_data)
}
}
