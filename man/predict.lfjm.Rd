\name{predict.lfjm}
\alias{predict.lfjm}
\title{Obtain Predictions from an LFJM Fit}
\description{
Extracts in-sample longitudinal fitted values, fitted functional scores, or
subject-specific random intercepts, and computes landmark-specific current
values and survival risk scores from a fitted longitudinal functional joint
model.
}
\usage{
\method{predict}{lfjm}(
  object,
  type = c(
    "longitudinal",
    "current_value",
    "risk",
    "random_effects",
    "scores"
  ),
  landmark = NULL,
  ...
)
}
\arguments{
  \item{object}{
  An object of class \code{"lfjm"}, typically returned by
  \code{\link{lfjm}}.
  }

  \item{type}{
  A character string specifying the requested prediction. The available
  choices are \code{"longitudinal"}, \code{"current_value"}, \code{"risk"},
  \code{"random_effects"}, and \code{"scores"}. Partial matching is allowed.
  }

  \item{landmark}{
  \code{NULL}, or a finite numeric scalar giving the landmark time.
  A landmark is required for \code{type = "current_value"} and
  \code{type = "risk"}. It is optional for
  \code{type = "random_effects"} and determines whether full-fit or
  history-based random intercepts are returned. It is ignored for
  \code{type = "longitudinal"} and \code{type = "scores"}.
  }

  \item{\dots}{
  Additional arguments. They are currently not used.
  }
}
\details{
The longitudinal functional joint model uses the conditional longitudinal
mean
\deqn{
\mu_{ij} =
\beta_0+\beta_{\mathrm{time}}t_{ij}+\beta_Q Q_{ij}
+\boldsymbol{\xi}_i^\mathsf{T}\boldsymbol{\theta}_B
+\boldsymbol{\zeta}_{ij}^\mathsf{T}\boldsymbol{\theta}_U
+u_i,
}{
mu_ij = beta_0 + beta_time t_ij + beta_Q Q_ij
+ xi_i' theta_B + zeta_ij' theta_U + u_i.
}
Here, \eqn{\boldsymbol{\xi}_i} and
\eqn{\boldsymbol{\zeta}_{ij}} are the between- and within-level scores from
Stage 1, while \eqn{u_i} is the shared subject-specific random intercept
profiled during Stage-2 adjusted profile h-likelihood estimation.

\code{type = "longitudinal"} returns the conditional fitted values
\eqn{\widehat{\mu}_{ij}} for the stored longitudinal observations. These
values use the original fitted scores and the full-fit random-intercept modes
in \code{object$uhat}. They include both between- and within-level functional
contributions.

\code{type = "scores"} returns the fitted Stage-1 score matrices without
modification. The between-level matrix has one row per subject, whereas the
within-level matrix has one row per longitudinal visit.

When \code{type = "random_effects"} and \code{landmark = NULL}, the method
returns \code{object$uhat}. These are the random-intercept modes obtained from
the full joint-model fit and therefore use all information entering that fit.

\subsection{Landmark predictions}{
For landmark time \eqn{L}, the method retains stored longitudinal and
functional observations satisfying \eqn{t_{ij}\leq L}. It then projects these
truncated functional histories onto the fitted Stage-1 mean surface and
eigensystems. Thus, the between-subject scores used for landmark prediction
are re-estimated from the functional history available by \eqn{L}, rather than
copied from the full-history fit.

Let \eqn{n_i(L)} be the number of longitudinal measurements available for
subject \eqn{i} by \eqn{L}, and let
\eqn{\widehat{\mu}^{\mathrm{fix}}_{ij,L}} denote the longitudinal predictor
formed from the scalar terms and the reprojected between- and within-level
scores, excluding the random intercept. The history-based random-intercept
predictor is
\deqn{
\widehat{u}_i(L) =
\frac{
\sum_{j:t_{ij}\leq L}
\{Y_{ij}-\widehat{\mu}^{\mathrm{fix}}_{ij,L}\}
}{
n_i(L)+\widehat{\sigma}_e^2/\widehat{\sigma}_u^2
}.
}{
u_hat_i(L) =
sum_{j:t_ij <= L}\{Y_ij - mu_fixed_hat_ij,L\} /
\{n_i(L) + sigma_e_hat^2 / sigma_u_hat^2\}.
}
This is the Gaussian random-intercept shrinkage predictor based on the
longitudinal history available by the landmark. It is returned when
\code{type = "random_effects"} and a landmark is supplied.

The shared current value at the landmark is
\deqn{
\widehat{m}_i(L) =
\widehat{\beta}_0
+\widehat{\beta}_{\mathrm{time}}L
+\widehat{\beta}_Q Q_i
+\widehat{\boldsymbol{\xi}}_i(L)^\mathsf{T}
 \widehat{\boldsymbol{\theta}}_B
+\widehat{u}_i(L).
}{
m_hat_i(L) = beta_0_hat + beta_time_hat L + beta_Q_hat Q_i
+ xi_i_hat(L)' theta_B_hat + u_i_hat(L).
}
Consistent with the joint-model specification, the transient within-visit
term does not enter the smooth current value used by the survival submodel.
This quantity is returned by \code{type = "current_value"}.

For \code{type = "risk"}, the returned value is the survival linear predictor
\deqn{
\widehat{r}_i(L) =
\widehat{\gamma}_0+\widehat{\gamma}_Q Q_i
+\widehat{\alpha}\widehat{m}_i(L).
}{
r_hat_i(L) = gamma_0_hat + gamma_Q_hat Q_i
+ alpha_hat m_i_hat(L).
}
It is a relative log-risk score. It does not include the Weibull baseline
log-hazard and is not a hazard, event probability, or survival probability.

Landmark results are produced for the subject identifiers stored in the
survival data. A subject with no longitudinal observation by \eqn{L} receives
\code{NA}. An error is raised if the requested landmark precedes all stored
longitudinal observations, so that no history is available for any subject.
The function does not restrict the returned subjects to those still at risk
at \eqn{L}.
}
}
\value{
The returned object depends on \code{type}:
\describe{
  \item{\code{"longitudinal"}}{
  An unnamed numeric vector with one conditional fitted value for each row of
  \code{object$longitudinal_data}.
  }
  \item{\code{"current_value"}}{
  A named numeric vector of landmark-specific current values, with names given
  by the subject identifiers in \code{object$survival_data}.
  }
  \item{\code{"risk"}}{
  A named numeric vector of landmark-specific relative log-risk scores, with
  names given by the subject identifiers in \code{object$survival_data}.
  }
  \item{\code{"random_effects"}}{
  A named numeric vector of subject-specific random intercepts. With
  \code{landmark = NULL}, these are the full-fit APHL modes. With a supplied
  landmark, they are the history-based shrinkage predictors
  \eqn{\widehat{u}_i(L)}.
  }
  \item{\code{"scores"}}{
  A list with components \code{between} and \code{within}. These contain,
  respectively, \code{object$lfpca$scoreB_subject} and
  \code{object$lfpca$scoreU_visit}.
  }
}
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

longitudinal_predictions <- predict(
  fit,
  type = "longitudinal"
)

full_random_intercepts <- predict(
  fit,
  type = "random_effects"
)

current_values <- predict(
  fit,
  type = "current_value",
  landmark = 2
)

risk_scores <- predict(
  fit,
  type = "risk",
  landmark = 2
)

landmark_random_intercepts <- predict(
  fit,
  type = "random_effects",
  landmark = 2
)

functional_scores <- predict(
  fit,
  type = "scores"
)

head(longitudinal_predictions)
head(current_values)
head(risk_scores)
head(landmark_random_intercepts)
dim(functional_scores$between)
dim(functional_scores$within)
}
}
