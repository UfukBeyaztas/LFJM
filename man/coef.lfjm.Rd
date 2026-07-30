\name{coef.lfjm}
\alias{coef.lfjm}
\title{Extract Coefficients from an LFJM Fit}
\description{
Extracts scalar parameters, score-space coefficients, and reconstructed
functional coefficient estimates from a fitted longitudinal functional joint
model.
}
\usage{
\method{coef}{lfjm}(
  object,
  type = c("all", "scalar", "between", "within", "survival", "scores"),
  ...
)
}
\arguments{
  \item{object}{
    An object of class \code{"lfjm"}, typically returned by
    \code{\link{lfjm}}.
  }
  \item{type}{
    A character string specifying the coefficients to return. The available
    choices are \code{"all"}, \code{"scalar"}, \code{"between"},
    \code{"within"}, \code{"survival"}, and \code{"scores"}. Partial
    matching is allowed.
  }
  \item{\dots}{
    Additional arguments. They are currently not used.
  }
}
\details{
The fitted model separates the centered longitudinal functional predictor into
a subject-level component, \eqn{B_i(s)}, and a visit-level component,
\eqn{U_{ij}(s)}. In the random-intercept implementation, the longitudinal
submodel is
\deqn{
Y_{ij} =
\beta_0 + \beta_{\mathrm{time}}t_{ij} + \beta_Q Q_i
+ \langle B_i,\beta_B\rangle
+ \langle U_{ij},\beta_U\rangle
+ u_i + \varepsilon_{ij}.
}{
Y_ij = beta_0 + beta_time t_ij + beta_Q Q_i
+ <B_i,beta_B> + <U_ij,beta_U> + u_i + epsilon_ij.
}
Here, \eqn{\beta_B} describes the effect of persistent between-subject
functional variation, whereas \eqn{\beta_U} describes the effect of
visit-specific within-subject variation.

The current value entering the survival model excludes the visit-specific
functional deviation:
\deqn{
m_i(t) =
\beta_0 + \beta_{\mathrm{time}}t + \beta_Q Q_i
+ \langle B_i,\beta_B\rangle + u_i.
}{
m_i(t) = beta_0 + beta_time t + beta_Q Q_i
+ <B_i,beta_B> + u_i.
}
The Weibull proportional-hazards submodel is
\deqn{
h_i(t) =
\tau t^{\tau-1}
\exp\{\gamma_0+\gamma_Q Q_i+\alpha m_i(t)\}.
}{
h_i(t) = tau t^(tau-1)
exp\{gamma_0 + gamma_Q Q_i + alpha m_i(t)\}.
}
Consequently, the functional effect induced in the survival linear predictor
is \eqn{\alpha\beta_B(s)}, rather than a separate direct functional effect.

At Stage 1, longitudinal functional principal component analysis supplies
between- and within-level eigenfunctions and predicted scores. At Stage 2,
the adjusted profile h-likelihood estimates the retained score coefficients
\eqn{\theta_{B,k}} and \eqn{\theta_{U,l}}. The reported functional
coefficients are reconstructed on the fitted functional grid as
\deqn{
\widehat{\beta}_B(s)
= \sum_{k=1}^{K_B}\widehat{\theta}_{B,k}
\widehat{\phi}^{B}_k(s),
\qquad
\widehat{\beta}_U(s)
= \sum_{l=1}^{K_U}\widehat{\theta}_{U,l}
\widehat{\phi}^{U}_l(s),
}{
beta_B_hat(s) = sum_k theta_B_hat[k] phi_B_hat[k](s),
beta_U_hat(s) = sum_l theta_U_hat[l] phi_U_hat[l](s).
}
The survival functional effect returned by \code{type = "survival"} is
\eqn{\widehat{\alpha}\widehat{\beta}_B(s)}.

The possible values of \code{type} have the following meanings:
\describe{
  \item{\code{"scalar"}}{
    Returns the finite-dimensional parameters
    \code{beta0}, \code{beta_time}, \code{beta_Q}, \code{sigma_e},
    \code{sigma_u}, \code{gamma0}, \code{gamma_Q}, \code{alpha}, and
    \code{Tau}. The quantities \code{sigma_e} and \code{sigma_u} are
    standard deviations.
  }
  \item{\code{"between"}}{
    Returns \eqn{\widehat{\beta}_B(s)} evaluated on the fitted functional
    grid.
  }
  \item{\code{"within"}}{
    Returns \eqn{\widehat{\beta}_U(s)} evaluated on the fitted functional
    grid.
  }
  \item{\code{"survival"}}{
    Returns
    \eqn{\widehat{\alpha}\widehat{\beta}_B(s)} evaluated on the fitted
    functional grid.
  }
  \item{\code{"scores"}}{
    Returns the retained between- and within-level score coefficients.
  }
  \item{\code{"all"}}{
    Returns all of the preceding coefficient groups in a single list.
  }
}

As in the fitted model, these estimates are conditional on the Stage-1 mean
surface, eigensystems, and predicted functional scores.
}
\value{
The returned object depends on \code{type}:
\describe{
  \item{\code{"scalar"}}{
    A named numeric vector with elements \code{beta0}, \code{beta_time},
    \code{beta_Q}, \code{sigma_e}, \code{sigma_u}, \code{gamma0},
    \code{gamma_Q}, \code{alpha}, and \code{Tau}.
  }
  \item{\code{"between"}}{
    A numeric vector containing the estimated between-level coefficient
    function on the fitted grid.
  }
  \item{\code{"within"}}{
    A numeric vector containing the estimated within-level coefficient
    function on the fitted grid.
  }
  \item{\code{"survival"}}{
    A numeric vector containing the estimated induced survival functional
    effect on the fitted grid.
  }
  \item{\code{"scores"}}{
    A list with components \code{between} and \code{within}, containing the
    estimated score-space coefficient vectors.
  }
  \item{\code{"all"}}{
    A list with the following components:
    \describe{
      \item{\code{scalar}}{
        The named vector of scalar parameter estimates.
      }
      \item{\code{score_coefficients}}{
        A list containing the between- and within-level score coefficients.
      }
      \item{\code{beta_between}}{
        The estimated between-level functional coefficient.
      }
      \item{\code{beta_within}}{
        The estimated within-level functional coefficient.
      }
      \item{\code{survival_between}}{
        The induced survival functional effect
        \eqn{\widehat{\alpha}\widehat{\beta}_B(s)}.
      }
    }
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

coef(fit)
coef(fit, type = "scalar")
coef(fit, type = "between")
coef(fit, type = "within")
coef(fit, type = "survival")
coef(fit, type = "scores")
}
}
