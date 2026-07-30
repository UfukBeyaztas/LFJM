\name{simulate_lfjm}
\alias{simulate_lfjm}
\title{Simulate Longitudinal Functional Joint-Model Data}
\description{
Generates data from one of two fixed longitudinal functional joint-model
settings: a longitudinal setting with visit-specific functional variation or
a baseline setting in which one functional curve is repeated across visits.
}
\usage{
simulate_lfjm(
  setting = c("longitudinal", "baseline"),
  number_subjects = 400,
  number_visits = 8,
  grid = seq(0, 1, length.out = 101),
  visit_grid = seq(0, 5, length.out = number_visits),
  seed = NULL
)
}
\arguments{
  \item{setting}{
  A character string specifying the simulation setting. Use
  \code{"longitudinal"} for a repeatedly observed functional predictor with
  between- and within-subject variation, or \code{"baseline"} for a single
  subject-specific curve copied across visits. Partial matching is allowed.
  }

  \item{number_subjects}{
  A positive integer giving the number of subjects.
  }

  \item{number_visits}{
  A positive integer giving the number of planned visits per subject before
  event-time and censoring restrictions are applied.
  }

  \item{grid}{
  A finite, strictly increasing numeric vector giving the common functional
  observation grid. At least two grid points are required.
  }

  \item{visit_grid}{
  A numeric vector giving the planned visit times. Its length must equal
  \code{number_visits}. The default gives equally spaced visits from zero to
  five.
  }

  \item{seed}{
  An optional value supplied to \code{\link{set.seed}} before data
  generation. If \code{NULL}, the current random-number-generator state is
  used without resetting it.
  }
}
\details{
The function implements the two data-generating settings used to illustrate
the longitudinal functional joint model. Let \eqn{s} denote the functional
domain and let \eqn{t_{ij}} be the planned visit time for subject \eqn{i} at
visit \eqn{j}. In the longitudinal setting, the observed functional predictor
is generated as
\deqn{
X_{ij}(s)
=
\eta(s,t_{ij})+B_i(s)+U_{ij}(s)+e_{ij}(s),
}{
X_ij(s) = eta(s,t_ij) + B_i(s) + U_ij(s) + e_ij(s),
}
where
\deqn{
\eta(s,t)
=
0.25\sin(\pi s)
+0.10\widetilde{t}\cos(2\pi s)
}{
eta(s,t) = 0.25 sin(pi s) + 0.10 t_tilde cos(2 pi s)
}
and \eqn{\widetilde{t}} is visit time standardized over the planned visit
grid. The latent components have finite Karhunen--Loeve representations
\deqn{
B_i(s)=\sum_{k=1}^{3}\xi_{ik}\phi_k^B(s),
\qquad
U_{ij}(s)=\sum_{l=1}^{2}\zeta_{ijl}\phi_l^U(s).
}{
B_i(s) = sum_k xi_ik phi_Bk(s),
U_ij(s) = sum_l zeta_ijl phi_Ul(s).
}
The trigonometric basis functions are orthonormalized using trapezoidal
integration weights on \code{grid}. The mutually independent scores satisfy
\eqn{\xi_{ik}\sim N(0,\lambda_k^B)} and
\eqn{\zeta_{ijl}\sim N(0,\lambda_l^U)}. The between-level eigenvalues are
\eqn{(3.0,1.5,0.6)} in both settings. In the longitudinal setting, the
within-level eigenvalues are \eqn{(0.5,0.25)} and the grid-error standard
deviation is \eqn{0.015}.

The scalar longitudinal response is generated from
\deqn{
Y_{ij}
=
\beta_0+\beta_{\mathrm{time}}t_{ij}+\beta_Q Q_i
+\boldsymbol{\xi}_i^\mathsf{T}\boldsymbol{\theta}_B
+\boldsymbol{\zeta}_{ij}^\mathsf{T}\boldsymbol{\theta}_U
+u_i+\varepsilon_{ij},
}{
Y_ij = beta_0 + beta_time t_ij + beta_Q Q_i
+ xi_i' theta_B + zeta_ij' theta_U + u_i + epsilon_ij,
}
where \eqn{Q_i\sim\operatorname{Bernoulli}(0.5)},
\eqn{u_i\sim N(0,\sigma_u^2)}, and
\eqn{\varepsilon_{ij}\sim N(0,\sigma_\varepsilon^2)}. The coefficient
vectors are
\eqn{\boldsymbol{\theta}_B=(0.90,-0.55,0)^\mathsf{T}} and, in the
longitudinal setting,
\eqn{\boldsymbol{\theta}_U=(0.45,-0.30)^\mathsf{T}}.

The survival submodel is linked to the smooth subject-specific current value
\deqn{
m_i(t)
=
\beta_0+\beta_{\mathrm{time}}t+\beta_Q Q_i
+\boldsymbol{\xi}_i^\mathsf{T}\boldsymbol{\theta}_B+u_i
}{
m_i(t) = beta_0 + beta_time t + beta_Q Q_i
+ xi_i' theta_B + u_i
}
through the Weibull proportional-hazards model
\deqn{
h_i(t)
=
\tau t^{\tau-1}
\exp\{\gamma_0+\gamma_Q Q_i+\alpha m_i(t)\}.
}{
h_i(t) = tau t^(tau-1)
exp\{gamma_0 + gamma_Q Q_i + alpha m_i(t)\}.
}
The within-visit functional deviation does not enter the hazard directly. Its
effect is confined to the longitudinal response, whereas the functional
effect induced in the log hazard is
\eqn{\alpha\beta_B(s)}.

The common scalar parameter values are
\deqn{
(\beta_0,\beta_{\mathrm{time}},\beta_Q,\sigma_u,
\sigma_\varepsilon,\tau,\gamma_0,\gamma_Q,\alpha)
=
(1.00,0.28,0.45,0.40,0.30,1.20,-4.20,-0.45,0.90).
}{
(beta_0, beta_time, beta_Q, sigma_u, sigma_epsilon, tau,
gamma_0, gamma_Q, alpha)
= (1.00, 0.28, 0.45, 0.40, 0.30, 1.20, -4.20, -0.45, 0.90).
}
Event times are obtained by numerical inversion of the cumulative hazard.
The censoring time is the minimum of an exponential random variable with
rate \eqn{0.05} and administrative censoring at time \eqn{8}. Only planned
visits satisfying \eqn{t_{ij}\leq T_i}, where \eqn{T_i} is the observed event
or censoring time, are retained. Consequently, the returned number of
longitudinal records can be smaller than
\code{number_subjects * number_visits}.

For \code{setting = "baseline"}, the within-level coefficient is zero and
the first retained functional curve for each subject is copied exactly to all
of that subject's retained visits. The core generator uses negligible
within-level eigenvalues \eqn{(10^{-8},10^{-8})} and grid-error standard
deviation \eqn{10^{-4}} before the first curve is copied. The returned true
within-subject curves and within-level coefficient function are set exactly
to zero. This setting represents the degenerate case in which the functional
predictor is time invariant.
}
\value{
An object of class \code{c("lfjm_simulation", "list")} with the following
components:
\describe{
  \item{\code{data.long}}{
  A data frame containing the retained longitudinal records, ordered by
  subject and visit time. It contains \code{ID}, \code{obstime},
  \code{row_id}, \code{Q1}, \code{time}, \code{event}, \code{Y}, and
  \code{true_mu_BU}. The observed functional predictor is stored in columns
  \code{func.X.1}, \code{func.X.2}, and so on. The variables \code{time} and
  \code{event} repeat the subject-level survival outcome on each retained
  visit. The variable \code{true_mu_BU} is the true conditional
  longitudinal mean, including the random intercept.
  }

  \item{\code{data.surv}}{
  A data frame containing one row per subject and the variables \code{ID},
  observed follow-up \code{time}, event indicator \code{event}, scalar
  covariate \code{Q1}, and true random intercept \code{true_u}.
  }

  \item{\code{grid}, \code{sgrid}}{
  Identical numeric vectors containing the functional grid.
  }

  \item{\code{true_B_visit}}{
  A numeric matrix containing the true between-subject functional component
  for each retained row of \code{data.long}. Its columns correspond to
  \code{grid}.
  }

  \item{\code{true_U_visit}}{
  A numeric matrix containing the true visit-specific within-subject
  functional component for each retained row of \code{data.long}. It is an
  all-zero matrix in the baseline setting.
  }

  \item{\code{truth}}{
  A list containing the generating values \code{beta0}, \code{beta_time},
  \code{beta_Q}, \code{sigma_u}, \code{sigma_eps}, \code{Tau},
  \code{gamma0}, \code{gamma_Q}, and \code{alpha}; score coefficients
  \code{theta_B} and \code{theta_U}; grid-evaluated coefficient functions
  \code{beta_B_true} and \code{beta_U_true}; the induced survival functional
  effect \code{surv_B_current_true}; basis matrices \code{phi0} and
  \code{phiU}; and the subject-level between scores \code{xi}.
  }

  \item{\code{case_type}}{
  A character string equal to \code{"longitudinal"} or \code{"baseline"}.
  }
}
}
\examples{
\dontrun{
simulated_data <- simulate_lfjm(
  setting = "longitudinal",
  number_subjects = 100,
  number_visits = 6,
  seed = 2026
)

print(simulated_data)
names(simulated_data$truth)

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
summary(fit)

baseline_data <- simulate_lfjm(
  setting = "baseline",
  number_subjects = 100,
  number_visits = 6,
  seed = 2026
)

print(baseline_data)
}
}
