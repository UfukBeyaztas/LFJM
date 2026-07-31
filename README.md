# LFJM: Longitudinal Functional Joint Models

[![R](https://img.shields.io/badge/R-%3E%3D%203.5.0-276DC3.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/UfukBeyaztas/LFJM)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-green.svg)](https://www.gnu.org/licenses/gpl-3.0)

`LFJM` fits joint models for studies in which a functional predictor is
observed repeatedly together with a scalar longitudinal response and a
time-to-event outcome.

The package separates persistent between-subject functional variation from
visit-specific within-subject variation using longitudinal functional
principal component analysis (LFPCA). The resulting longitudinal-survival
model is estimated by adjusted profile h-likelihood (APHL) under a shared
current-value association.

The complete function reference is available in the
[LFJM 1.0.0 package manual](LFJM_1.0.0.pdf).

## Methodological overview

For subject \(i\) at visit \(j\), the repeatedly observed functional predictor
is decomposed as

$$
X_{ij}(s)=\eta(s,t_{ij})+B_i(s)+U_{ij}(s)+e_{ij}(s),
$$

where \(B_i(s)\) is the persistent between-subject component and
\(U_{ij}(s)\) is the visit-specific within-subject component. The scalar
longitudinal submodel is

$$
Y_{ij}
=\beta_0+\beta_{\mathrm{time}}t_{ij}+\beta_QQ_i
+\langle B_i,\beta_B\rangle
+\langle U_{ij},\beta_U\rangle
+u_i+\varepsilon_{ij}.
$$

The survival process is linked to the smooth current value

$$
m_i(t)
=\beta_0+\beta_{\mathrm{time}}t+\beta_QQ_i
+\langle B_i,\beta_B\rangle+u_i
$$

through a Weibull proportional-hazards model,

$$
h_i(t)=\tau t^{\tau-1}
\exp\{\gamma_0+\gamma_QQ_i+\alpha m_i(t)\}.
$$

The visit-specific component enters the longitudinal response through
\(\beta_U\), but it does not enter the smooth current value used by the
survival submodel. The induced functional effect in the survival linear
predictor is therefore \(\alpha\beta_B(s)\).

Estimation has two stages:

1. Weighted-grid LFPCA estimates the mean surface, between- and within-level
   eigensystems, and functional scores.
2. APHL estimates the longitudinal and survival parameters while profiling
   the subject-specific random intercepts.

Stage 2 treats the estimated Stage-1 quantities as fixed.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("UfukBeyaztas/LFJM")
```

Then load the package:

```r
library(LFJM)
```

## Data structure

The main fitting function requires four aligned objects:

| Object | Required structure |
|---|---|
| `longitudinal_data` | One row per retained visit. By default, the required columns are `ID`, `obstime`, `Y`, and `Q1`. |
| `survival_data` | Exactly one row per subject. By default, the required columns are `ID`, `time`, `event`, and `Q1`. |
| `functional_data` | A numeric matrix with one row per row of `longitudinal_data` and one column per functional grid point. |
| `grid` | A finite, strictly increasing numeric vector whose length equals `ncol(functional_data)`. |

Alternative column names can be supplied through the corresponding arguments
of `lfjm()`.

## Quick start

The package includes a generator for longitudinal and baseline functional
joint-model settings.

```r
library(LFJM)

simulated_data <- simulate_lfjm(
  setting = "longitudinal",
  number_subjects = 100,
  number_visits = 6,
  seed = 2026
)

print(simulated_data)
```

Extract the functional matrix from the simulated longitudinal data:

```r
functional_columns <- paste0(
  "func.X.",
  seq_along(simulated_data$grid)
)

functional_data <- as.matrix(
  simulated_data$data.long[
    ,
    functional_columns,
    drop = FALSE
  ]
)
```

Specify the estimation controls and fit the model:

```r
control <- lfjm_control(
  pve_between = 0.95,
  pve_within = 0.95,
  presmooth = TRUE,
  smoothing_df = 12,
  quadrature_points = 32,
  max_iterations = 500,
  multistart = TRUE,
  standard_errors = TRUE
)

fit <- lfjm(
  longitudinal_data = simulated_data$data.long,
  survival_data = simulated_data$data.surv,
  functional_data = functional_data,
  grid = simulated_data$grid,
  control = control
)

print(fit)
summary(fit)
```

When desired, fixed component counts can be supplied through
`components_between` and `components_within` in `lfjm_control()`. Otherwise,
the component counts are selected using the specified proportions of variance
explained.

## Extracting model results

Scalar parameters, score coefficients, and reconstructed functional
coefficients are available through `coef()`:

```r
scalar_parameters <- coef(fit, type = "scalar")
score_coefficients <- coef(fit, type = "scores")

beta_between <- coef(fit, type = "between")
beta_within <- coef(fit, type = "within")
survival_functional_effect <- coef(fit, type = "survival")

all_coefficients <- coef(fit, type = "all")
```

Additional fitted quantities can be obtained with standard R generics:

```r
longitudinal_fitted_values <- fitted(fit)

conditional_residuals <- residuals(
  fit,
  type = "response"
)

marginal_residuals <- residuals(
  fit,
  type = "marginal"
)

adjusted_log_likelihood <- logLik(fit)
conditional_covariance <- vcov(fit)
```

`vcov(fit)` returns the conditional Stage-2 covariance matrix on the
optimization scale. In particular, the matrix uses `log_sigma_e`,
`log_sigma_u`, and `log_Tau`; see the package manual for the natural-scale
delta-method transformation.

## Prediction

The fitted model supports longitudinal fitted values, functional scores,
random-intercept predictions, and landmark predictions:

```r
longitudinal_predictions <- predict(
  fit,
  type = "longitudinal"
)

functional_scores <- predict(
  fit,
  type = "scores"
)

full_fit_random_intercepts <- predict(
  fit,
  type = "random_effects"
)

current_values_at_2 <- predict(
  fit,
  type = "current_value",
  landmark = 2
)

risk_scores_at_2 <- predict(
  fit,
  type = "risk",
  landmark = 2
)

history_based_random_intercepts <- predict(
  fit,
  type = "random_effects",
  landmark = 2
)
```

The landmark risk output is a relative log-risk score,
\(\gamma_0+\gamma_QQ_i+\alpha m_i(L)\). It is not an event probability,
survival probability, or hazard.

## Visualization

Plot the reconstructed longitudinal and induced survival functional
coefficients:

```r
plot(fit, type = "coefficients")
```

Plot the retained between- and within-level eigensystems:

```r
plot(
  fit,
  type = "eigenfunctions",
  level = "between"
)

plot(
  fit,
  type = "eigenfunctions",
  level = "within"
)

plot(
  fit,
  type = "eigenvalues",
  level = "between"
)

plot(
  fit,
  type = "eigenvalues",
  level = "within"
)
```

## Main functions and methods

| Function or method | Purpose |
|---|---|
| `lfjm()` | Fit the two-stage longitudinal functional joint model. |
| `lfjm_control()` | Set LFPCA, quadrature, optimization, and standard-error controls. |
| `lfjm_lfpca()` | Fit the Stage-1 longitudinal functional principal component decomposition directly. |
| `simulate_lfjm()` | Generate data under longitudinal or baseline functional-predictor settings. |
| `coef()` | Extract scalar, score-space, between-level, within-level, or induced survival coefficients. |
| `summary()` | Construct the conditional Stage-2 coefficient table and fit summary. |
| `predict()` | Obtain longitudinal, score, random-effect, current-value, or landmark risk predictions. |
| `fitted()` | Extract conditional longitudinal fitted values. |
| `residuals()` | Extract conditional response or marginal longitudinal residuals. |
| `vcov()` | Extract the conditional Stage-2 inverse-Hessian covariance matrix. |
| `logLik()` | Extract the maximized adjusted profile h-likelihood. |
| `plot()` | Plot coefficient functions, eigenfunctions, or retained eigenvalues. |

For detailed arguments, returned components, equations, and interpretation,
consult the [package manual](LFJM_1.0.0.pdf) or use R help:

```r
help(package = "LFJM")
?lfjm
?lfjm_control
?lfjm_lfpca
?predict.lfjm
?simulate_lfjm
```

## Current model scope

Version 1.0.0 implements:

- a repeatedly observed functional predictor on a common grid;
- a Gaussian scalar longitudinal outcome;
- separate between- and within-subject functional effects;
- a subject-specific random intercept;
- a Weibull proportional-hazards survival submodel;
- a shared current-value association;
- one scalar longitudinal covariate and one scalar survival covariate; and
- conditional Stage-2 standard errors based on the numerical APHL Hessian.

Uncertainty from the estimated Stage-1 mean surface, eigensystems, and scores
is not propagated into the reported Stage-2 covariance matrix.

## Citation

To obtain the package citation in R, use:

```r
citation("LFJM")
```

When using the methodology in scientific work, please also cite the
accompanying methodological paper.

## License

`LFJM` is released under the
[GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0).

## Author

Ufuk Beyaztas  
Email: [ufukbeyaztas@gmail.com](mailto:ufukbeyaztas@gmail.com)  
ORCID: [0000-0002-5208-4950](https://orcid.org/0000-0002-5208-4950)

Questions and bug reports can be submitted through the
[GitHub issue tracker](https://github.com/UfukBeyaztas/LFJM/issues).
