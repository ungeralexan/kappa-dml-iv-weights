# Implementation notes

This document describes the computational implementation used for the thesis.
It complements the methodological discussion in `draft_of_13.tex`; software
details that are not needed to understand the estimators belong here rather
than in the main thesis text.

## Status of the code

The project currently uses one shared R source file, `functions_all.R`. It is a
project-level function library, not a formally installed R package: the
repository does not contain a `DESCRIPTION` file or `NAMESPACE`. Each
application notebook sources the same root file so that the estimator
definitions, standard errors, outcome-weight constructors, and diagnostics do
not drift across applications.

The shared code is used by the Vietnam, Card, and Angrist–Evans analyses. It
implements the classical kappa estimators and connects them to the DML and
outcome-weight workflow used throughout the thesis.

## Why a shared implementation was developed

When the thesis implementation was designed, the available software did not
provide the complete workflow required for the empirical comparison in one R
environment.

- `causalweight` provided an R implementation of an inverse-instrument-
  propensity weighting estimator for the LATE, with bootstrap-based inference,
  but not the complete normalized and unnormalized kappa menu used in the
  thesis together with the required covariate-balancing specification.
- The authors' `kappalate` implementation provided the complete estimator menu,
  maximum-likelihood and covariate-balancing instrument propensity scores, and
  analytical standard errors, but it was written for Stata.
- The thesis additionally required closed-form outcome weights for every kappa
  estimator, numerical reconstruction checks, common weight diagnostics,
  covariate-balance plots, and complete-pipeline outcome-transformation reruns
  alongside the DML estimators in R.

The custom implementation therefore placed the point estimators, propensity-
score procedures, analytical standard errors, outcome weights, and diagnostics
in one reproducible R workflow.

The software landscape changed during the project. Since June 2026, the CRAN
package `drlate` has implemented the complete kappa estimator suite in R,
including maximum-likelihood and covariate-balancing instrument propensity
scores and joint M-estimation standard errors. The thesis code remains useful
because it is the implementation on which the reported analysis was developed
and because it includes the thesis-specific outcome-weight reconstruction,
diagnostic, learner-comparison, and transformation-rerun pipeline. A future
validation step can compare the classical point estimates and standard errors
against `drlate` in addition to the existing replication checks.

## Notation and variable conventions

The shared functions use the following arguments:

| Argument | Meaning |
|---|---|
| `Y` | observed outcome |
| `D` | binary treatment |
| `Z` | binary instrument |
| `X` | kappa covariate design matrix, including an intercept |
| `p` | fitted instrument propensity score, $P(Z=1\mid X)$ |
| `X_df` | covariate data frame without an explicit intercept, used by the 2SLS wrapper |

The LaTeX thesis distinguishes the population propensity score
$p(X_i)$ from a fitted score $\widehat p_i=\widehat p(X_i)$. The R
functions receive the fitted vector as `p`.

## Instrument propensity-score estimation

Two logit-based procedures are implemented.

### Maximum likelihood

`logit_mle()` fits the binary instrument by logistic maximum likelihood. The
analytical-standard-error workflow uses `fit_logit_alpha()`, which standardizes
non-intercept covariates before fitting to improve numerical conditioning.
This standardization is a reparameterization and does not change the fitted
propensity scores.

### Covariate balancing

`fit_cbps_alpha()` solves

$$
\frac{1}{N}\sum_{i=1}^N
\frac{Z_i-p_i}{p_i(1-p_i)}X_i=0.
$$

This moment condition is algebraically equivalent to equality of the
instrument-group inverse-probability-weighted covariate moments. The solver is
initialized at the logit maximum-likelihood estimate and uses Newton updates
with a backtracking line search. It records convergence and the largest
remaining absolute balancing moment. Fitted probabilities are clipped to
$[10^{-8},1-10^{-8}]$ for numerical stability.

The intercept in `X` is essential for the finite-sample equivalence among the
covariate-balanced normalized estimators.

## Kappa point estimators

The shared implementation computes the following reported specifications:

| Code name | Thesis notation | Propensity score |
|---|---|---|
| `tau_cb_u` | $\widehat\tau_u^{\mathrm{CB}}$ | covariate balancing |
| `tau_ml_u` | $\widehat\tau_u^{\mathrm{ML}}$ | logit maximum likelihood |
| `tau_ml_a10` | $\widehat\tau_{a,10}^{\mathrm{ML}}$ | logit maximum likelihood |
| `tau_ml_a` | $\widehat\tau_a^{\mathrm{ML}}$ | logit maximum likelihood |
| `tau_ml_t` | $\widehat\tau_t^{\mathrm{ML}}=\widehat\tau_{a,1}^{\mathrm{ML}}$ | logit maximum likelihood |
| `tau_ml_a0` | $\widehat\tau_{a,0}^{\mathrm{ML}}$ | logit maximum likelihood |

With covariate balancing and an intercept, the relevant normalized variants
are numerically equivalent. The application tables therefore report
$\widehat\tau_u^{\mathrm{CB}}$ once instead of printing duplicate CBPS rows.
The estimator $\widehat\tau_a^{\mathrm{CB}}$ is not part of this equivalence.

`kappa_point_estimates()` is the canonical point-estimation interface.
`kappa_analytic_se_all()` adds the corresponding analytical standard errors.

## Analytical standard errors

The kappa standard errors use stacked M-estimation systems. Each system includes
the propensity-score moments, the estimator-specific numerator and denominator
moments, and the final ratio or contrast. `sandwich_se_mest()` computes the
sandwich covariance matrix, using a numerical central-difference Jacobian.

The maximum-likelihood propensity-score block uses the logit score
$(Z_i-p_i)X_i$. The covariate-balancing block uses
$(Z_i-p_i)X_i/\{p_i(1-p_i)\}$. This construction propagates uncertainty from
the estimated instrument propensity score into the standard error of the LATE
estimator.

## Outcome weights and identity checks

`kappa_outcome_weights()` constructs a closed-form outcome-weight vector for
each reported kappa estimator. A valid vector must satisfy

$$
\widehat\tau=\sum_{i=1}^N\omega_iY_i
$$

to the project's numerical tolerance. `check_weight_identity()` performs this
reconstruction check. Candidate vectors that fail the check are not used in
weight-based diagnostics.

The diagnostics are computed from unrounded weights and include normalization
masses, effective support, maximum absolute weights, sign composition,
cancellation, and covariate balance. Presentation tables apply rounding only
after these quantities have been computed.

## DML integration

The shared library also contains helper functions for PLR-IV and Wald-AIPW
objects, outcome-weight extraction, reconstruction checks, and balance plots.
The application notebooks compare linear or logistic nuisance models with
honest forests, Ranger, and XGBoost variants. Expensive XGBoost tuning is
performed in separate notebooks and exported as `.rds` pipeline artifacts.

For PLR-IV, the outcome, instrument, and treatment nuisance functions are all
implemented as regression tasks. The parametric specification therefore uses
linear regression for all three functions. For Wald-AIPW, the
instrument-specific outcome functions are regression tasks, while the
instrument propensity and instrument-specific treatment functions are binary
classification tasks. Its parametric specification combines linear outcome
regressions with logistic instrument and treatment models. The Ranger and
XGBoost specifications follow the same regression/classification mapping.

The headline outcome-weight specifications are fitted with
`OutcomeWeights::dml_with_smoother()` and honest
`grf::regression_forest()` nuisance learners. This implementation uses
regression forests for every nuisance function, including functions with a
binary response. All DML specifications use seed 42 and five outer
cross-fitting folds.

The headline honest forests pass `tune.parameters = "all"` to each GRF fit
inside its outer training sample. GRF then uses its internal out-of-bag tuning
procedure to select the sample fraction, `mtry`, minimum node size, honesty
fraction, honesty pruning rule, split-balance parameter `alpha`, and imbalance
penalty. This is distinct from the explicit inner cross-validation used for
XGBoost.

The XGBoost specifications use fold-specific nested tuning through
`DoubleML` and `mlr3tuning`. Within each of the five outer training samples,
each nuisance learner is tuned separately by three-fold inner
cross-validation and 15 random-search evaluations. The search varies
`nrounds` from 50 to 400, `max_depth` from 2 to 5, `eta` from 0.02 to 0.20,
and `min_child_weight` from 1 to 8. PLR-IV uses RMSE for all three nuisance
functions. Wald-AIPW uses RMSE for its outcome regression and Bernoulli
log-loss for its instrument and treatment classification models.
`tune_on_folds = TRUE` ensures that the outer evaluation fold is excluded
from both nuisance fitting and hyperparameter selection, and it implies that
the selected hyperparameters may differ across outer folds.

The Vietnam and Card tuning notebooks default to loading their verified
exports with `RERUN_TUNING = FALSE`. The current Child tuning notebooks have
no corresponding switch and rerun tuning when rendered; their stored exports
are consumed directly by the downstream Child analysis.

Complete-pipeline transformation checks rerun the relevant nuisance estimation,
cross-fitting, and tuning rules after the outcome has been shifted. They are
distinct from frozen-weight algebraic checks, which hold a fitted weight vector
fixed and evaluate the implied change directly.

## Reproducibility

The repository-level `README.md` records the validated R and package versions,
required data files, notebook execution order, stored tuning exports, and cache
policy. For a public replication repository, retain the canonical
`functions_all.R`, the application notebooks, the verified `.rds` exports, and
an `renv.lock` once it has been created. Do not distribute duplicate helper
files inside the application folders.

## External software references

- `kappalate`: Stata implementation accompanying Słoczyński, Uysal, and
  Wooldridge, with the kappa estimator menu, MLE/CBPS instrument propensity
  scores, and analytical inference.
- `causalweight`: R package containing inverse-probability-weighted and related
  causal estimators, historically used as the closest R reference for
  $\widehat\tau_u$.
- `drlate` (version 0.3.1, released 24 June 2026): R package implementing the
  complete Słoczyński--Uysal--Wooldridge kappa menu and joint M-estimation
  inference.
- `OutcomeWeights`: source of the general outcome-weight extraction framework
  used for the DML comparison.
