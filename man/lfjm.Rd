\name{lfjm}
\alias{lfjm}
\title{Fit a Longitudinal Functional Joint Model}
\description{
Fits a two-stage longitudinal functional joint model for a repeatedly observed
functional predictor, a scalar longitudinal outcome, and a time-to-event
outcome. Longitudinal functional principal component analysis separates
between- and within-subject functional variation, and the joint model is
estimated by adjusted profile h-likelihood.
}
\usage{
lfjm(
  longitudinal_data,
  survival_data,
  functional_data,
  grid,
  id = "ID",
  visit_time = "obstime",
  response = "Y",
  longitudinal_covariate = "Q1",
  survival_time = "time",
  event = "event",
  survival_covariate = "Q1",
  control = lfjm_control()
)
}
\arguments{
  \item{longitudinal_data}{
    A data frame containing one row for each longitudinal visit. It must
    contain the subject identifier, visit time, scalar longitudinal response,
    and scalar longitudinal covariate specified by \code{id},
    \code{visit_time}, \code{response}, and
    \code{longitudinal_covariate}, respectively.
  }
  \item{survival_data}{
    A data frame containing exactly one row for each subject. It must contain
    the subject identifier, observed event or censoring time, event indicator,
    and scalar survival covariate specified by \code{id},
    \code{survival_time}, \code{event}, and
    \code{survival_covariate}, respectively.
  }
  \item{functional_data}{
    A numeric matrix containing the repeatedly observed functional predictor.
    Each row must correspond to the same row of \code{longitudinal_data}, and
    each column must correspond to one point in \code{grid}.
  }
  \item{grid}{
    A finite, strictly increasing numeric vector giving the common functional
    observation grid. Its length must equal the number of columns of
    \code{functional_data}.
  }
  \item{id}{
    A character string naming the subject-identifier column in both data
    frames.
  }
  \item{visit_time}{
    A character string naming the visit-time column in
    \code{longitudinal_data}.
  }
  \item{response}{
    A character string naming the scalar longitudinal-response column in
    \code{longitudinal_data}.
  }
  \item{longitudinal_covariate}{
    A character string naming the scalar covariate included in the
    longitudinal submodel.
  }
  \item{survival_time}{
    A character string naming the observed event or censoring-time column in
    \code{survival_data}. All values must be positive and finite.
  }
  \item{event}{
    A character string naming the event-indicator column in
    \code{survival_data}. The column must contain only \code{0} for censoring
    and \code{1} for an observed event.
  }
  \item{survival_covariate}{
    A character string naming the scalar covariate entering the survival
    submodel directly. It may refer to the same baseline covariate used in the
    longitudinal submodel.
  }
  \item{control}{
    An object returned by \code{\link{lfjm_control}}, specifying LFPCA,
    quadrature, optimization, and standard-error settings.
  }
}
\details{
The method is designed for a functional predictor observed repeatedly for each
subject. After removal of the estimated mean surface, the predictor is written
as
\deqn{
X_{ij}(s) = \eta(s,t_{ij}) + B_i(s) + U_{ij}(s) + e_{ij}(s),
}{
X_ij(s) = eta(s,t_ij) + B_i(s) + U_ij(s) + e_ij(s),
}
where \eqn{B_i(s)} represents persistent between-subject variation and
\eqn{U_{ij}(s)} represents visit-specific within-subject variation.

Estimation proceeds in two stages. At Stage 1,
\code{\link{lfjm_lfpca}} estimates the mean surface, the between- and
within-level eigensystems, and their scores. With retained scores
\eqn{\widehat{\boldsymbol{\xi}}_i} and
\eqn{\widehat{\boldsymbol{\zeta}}_{ij}}, Stage 2 fits the working
longitudinal submodel
\deqn{
Y_{ij} =
\beta_0 + \beta_{\mathrm{time}}t_{ij} + \beta_Q Q_i
+ \widehat{\boldsymbol{\xi}}_i^\mathsf{T}\boldsymbol{\theta}_B
+ \widehat{\boldsymbol{\zeta}}_{ij}^\mathsf{T}\boldsymbol{\theta}_U
+ u_i + \varepsilon_{ij},
}{
Y_ij = beta_0 + beta_time t_ij + beta_Q Q_i
+ xi_hat_i' theta_B + zeta_hat_ij' theta_U + u_i + epsilon_ij.
}
The subject-level current value is
\deqn{
m_i(t) =
\beta_0 + \beta_{\mathrm{time}}t + \beta_Q Q_i
+ \widehat{\boldsymbol{\xi}}_i^\mathsf{T}\boldsymbol{\theta}_B
+ u_i.
}{
m_i(t) = beta_0 + beta_time t + beta_Q Q_i
+ xi_hat_i' theta_B + u_i.
}
The visit-specific functional deviation is excluded from the current value
because it represents transient state variation. The survival submodel is
\deqn{
h_i(t) =
\tau t^{\tau-1}
\exp\{\gamma_0+\gamma_Q Q_i+\alpha m_i(t)\},
\qquad \tau>0.
}{
h_i(t) = tau t^(tau-1)
exp\{gamma_0 + gamma_Q Q_i + alpha m_i(t)\}, tau > 0.
}
Thus, the induced functional effect in the survival linear predictor is
\eqn{\alpha\beta_B(s)}.

The Stage-2 criterion is the adjusted profile h-likelihood
\deqn{
\widetilde{h}(\boldsymbol{\psi}) =
\sum_i
\left[
h_i\{\widehat{u}_i(\boldsymbol{\psi});\boldsymbol{\psi}\}
- \frac{1}{2}
\log\left\{\frac{H_i(\boldsymbol{\psi})}{2\pi}\right\}
\right],
}{
h_tilde(psi) = sum_i [h_i(u_hat_i(psi); psi)
- 0.5 log\{H_i(psi)/(2 pi)\}].
}
Here, \eqn{\widehat{u}_i} is the profiled random-intercept mode and
\eqn{H_i} is the negative second derivative of the subject-specific
h-likelihood with respect to that random intercept. The cumulative-hazard
integrals are evaluated by Gauss--Legendre quadrature after a transformation
that regularizes the Weibull integrand.

After estimation, the retained score coefficients are mapped back to the
functional grid:
\deqn{
\widehat{\beta}_B(s)
= \sum_{k=1}^{K_B}
\widehat{\theta}_{B,k}\widehat{\phi}^{B}_k(s),
\qquad
\widehat{\beta}_U(s)
= \sum_{l=1}^{K_U}
\widehat{\theta}_{U,l}\widehat{\phi}^{U}_l(s).
}{
beta_B_hat(s) = sum_k theta_B_hat[k] phi_B_hat[k](s),
beta_U_hat(s) = sum_l theta_U_hat[l] phi_U_hat[l](s).
}

Before estimation, \code{longitudinal_data} and the rows of
\code{functional_data} are jointly ordered by subject identifier and visit
time. The survival data are then matched to the ordered subject identifiers.
The subject identifiers appearing in the two data frames must agree exactly.
}
\value{
An object of class \code{"lfjm"}. It is a list containing the following
principal components:
\describe{
  \item{\code{method}}{
    A character string identifying the LFPCA--APHL estimator.
  }
  \item{\code{beta}}{
    The estimated longitudinal coefficient vector, including the intercept,
    visit-time effect, scalar-covariate effect, and retained between- and
    within-level score coefficients.
  }
  \item{\code{theta_B}, \code{theta_U}}{
    The retained between- and within-level score coefficient vectors.
  }
  \item{\code{beta_B}, \code{beta_U}}{
    The reconstructed between- and within-level functional coefficients
    evaluated on \code{grid}.
  }
  \item{\code{surv_B_current}}{
    The induced survival functional effect
    \eqn{\widehat{\alpha}\widehat{\beta}_B(s)} evaluated on \code{grid}.
  }
  \item{\code{se_beta_B}, \code{se_beta_U}, \code{se_surv_B}}{
    Pointwise conditional standard errors for the reconstructed functional
    coefficients and induced survival effect. These components are
    \code{NULL} if standard errors were not requested or could not be
    computed.
  }
  \item{\code{sigma_e}, \code{sigma_u}}{
    The estimated residual and random-intercept standard deviations.
  }
  \item{\code{gamma0}, \code{gamma_Q}, \code{alpha}, \code{Tau}}{
    The estimated survival intercept, direct scalar-covariate effect,
    current-value association parameter, and Weibull shape parameter.
  }
  \item{\code{uhat}}{
    The profiled subject-specific random-intercept modes.
  }
  \item{\code{risk}}{
    The fitted subject-specific survival linear predictors evaluated at the
    observed event or censoring times.
  }
  \item{\code{loglik}, \code{objective}, \code{npar}}{
    The maximized adjusted profile h-likelihood, minimized negative criterion,
    and number of fitted Stage-2 parameters.
  }
  \item{\code{parameter_vector}}{
    The complete optimization parameter vector. Positive scale parameters are
    represented on the log scale in this vector.
  }
  \item{\code{covariance_matrix}}{
    The conditional Stage-2 covariance matrix obtained from the numerical
    Hessian. It is \code{NULL} when unavailable.
  }
  \item{\code{convergence}, \code{convergence_message}}{
    The optimizer convergence code and optional convergence message.
  }
  \item{\code{lfpca}}{
    The fitted Stage-1 object of class \code{"lfjm_lfpca"}.
  }
  \item{\code{longitudinal_data}, \code{survival_data},
        \code{functional_data}, \code{grid}}{
    The internally ordered model data and functional grid.
  }
  \item{\code{call}, \code{control}, \code{variable_names}}{
    The matched call, control settings, and original variable-name mapping.
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

control <- lfjm_control(
  pve_between = 0.95,
  pve_within = 0.95,
  quadrature_points = 32,
  max_iterations = 500
)

fit <- lfjm(
  longitudinal_data = simulated_data$data.long,
  survival_data = simulated_data$data.surv,
  functional_data = functional_data,
  grid = simulated_data$grid,
  id = "ID",
  visit_time = "obstime",
  response = "Y",
  longitudinal_covariate = "Q1",
  survival_time = "time",
  event = "event",
  survival_covariate = "Q1",
  control = control
)

print(fit)
summary(fit)
coef(fit, type = "scalar")
fitted(fit)
predict(fit, type = "risk", landmark = 2)
}
}
