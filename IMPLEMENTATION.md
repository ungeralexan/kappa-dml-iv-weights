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
each reported kappa estimator from the observed instrument `Z`, treatment `D`,
and fitted instrument propensity scores `p`. Propensity-score estimation is
completed before this function is called; no outcome regression or additional
optimization is performed during kappa weight construction.

For the normalized Uysal estimator, the function computes

$$
s_1=\sum_i\frac{Z_i}{p_i},\qquad
s_0=\sum_i\frac{1-Z_i}{1-p_i},
$$

the normalized first-stage contrast

$$
\widehat\Delta_D^u
=
\frac{\sum_iD_iZ_i/p_i}{s_1}
-
\frac{\sum_iD_i(1-Z_i)/(1-p_i)}{s_0},
$$

and the final outcome weights

$$
\omega_i^u
=
\frac{
Z_i/(p_i s_1)
-
(1-Z_i)/\{(1-p_i)s_0\}
}{
\widehat\Delta_D^u
}.
$$

This is the formula derived in the thesis appendix. The other kappa vectors
are constructed as

$$
\omega_i^{a,10}
=
\frac{\kappa_{1i}}{\sum_j\kappa_{1j}}
-
\frac{\kappa_{0i}}{\sum_j\kappa_{0j}},
$$

while the three unnormalized estimators use

$$
\omega_i^a
=
\frac{(Z_i-p_i)/\{p_i(1-p_i)\}}{\sum_j\kappa_j},
\qquad
\omega_i^{a,1}
=
\frac{(Z_i-p_i)/\{p_i(1-p_i)\}}{\sum_j\kappa_{1j}},
\qquad
\omega_i^{a,0}
=
\frac{(Z_i-p_i)/\{p_i(1-p_i)\}}{\sum_j\kappa_{0j}}.
$$

In the R implementation these are returned as `w_a`, `w_a1`, and `w_a0`.

Every candidate vector must satisfy

$$
\widehat\tau=\sum_{i=1}^N\omega_iY_i
$$

`check_weight_identity()` performs this numerical reconstruction check. It
compares the unrounded inner product `sum(w * Y)` with the corresponding
unrounded point estimate by calling R's `all.equal()` without overriding its
default tolerance. The check is therefore scale-aware and uses R's default
numerical comparison rule (approximately the square root of machine precision,
rather than a fixed absolute threshold of exactly $10^{-8}$). This is the same
`all.equal()` convention reported in the application notebooks.

The reconstruction check is distinct from the normalization diagnostics. A
weight vector can reproduce its fitted estimate without having overall,
treated, and control masses equal to $0$, $1$, and $1$, respectively. Conversely,
a vector with approximately normalized masses is not treated as an exact
outcome-weight representation unless it also passes the reconstruction check.

### DML outcome-weight extraction

The DML outcome weights are obtained from the development version
`0.2.0.9000` of `OutcomeWeights`, using GitHub commit
`0c94f940b04c14d0247b46842af37752e306b79e`.

The headline PLR-IV and Wald-AIPW specifications use
`OutcomeWeights::dml_with_smoother()` with five-fold cross-fitting and honest
`grf` regression-forest smoothers. Their fitted observation-level vectors are
recovered with `OutcomeWeights::get_outcome_weights()`. The alternative
parametric, Ranger, and XGBoost specifications are fitted as
`DoubleMLPLIV` or `DoubleMLIIVM` objects. Where the learner and package method
support weight recovery, their stored nuisance fits are passed to the
corresponding `DoubleML` method of `get_outcome_weights()`.

The same `check_weight_identity()` helper is applied to weight vectors
extracted from `DoubleML` fits; no separate DoubleML-specific identity helper
is used. A fitted vector is classified as an exact outcome-weight
representation only when its inner product with the observed outcome passes
R's default `all.equal()` comparison against the fitted coefficient. Vectors
that fail this identity, including the affected XGBoost extractions in the
current analysis, may be retained as explicitly flagged learner-sensitivity
information but are excluded from the main normalization, concentration, and
covariate-balance comparison.

The diagnostics are computed from unrounded weights and include normalization
masses, effective support, maximum absolute weights, sign composition,
cancellation, and covariate balance. Presentation tables apply rounding only
after these quantities have been computed.

### Covariate-balance computation

All three application notebooks call the shared `make_love()` helper in
`functions_all.R`. The helper receives the treatment indicator, the
application-specific balance matrix, and a numerically verified outcome-weight
vector as explicit arguments. It first orients the signed outcome weights as
`w * (2 * D - 1)`, so treated coefficients retain their sign and control
coefficients are sign-reversed. `cobalt::love.plot()` then normalizes these
oriented weights separately within the observed treatment groups when
computing the adjusted covariate means.

The validated environment uses `cobalt` 4.6.2 with
`continuous = "std"`, `binary = "raw"`, `s.d.denom = "pooled"`,
`stats = "mean.diffs"`, and `abs = TRUE`. For a continuous covariate, the
adjusted and unadjusted mean differences are divided by the same unweighted
pooled within-group standard deviation,
`sqrt((var_T + var_C) / 2)`. For a binary covariate, `cobalt` detects the
two-valued column and reports the absolute raw difference between the two
signed weighted means without standardization. Consequently, the plotted
reference value of 0.1 represents 0.1 pooled standard deviations for a
continuous covariate but a raw signed-moment difference of 0.10 for a binary
covariate; it is a visual guide rather than a common pass/fail scale.

The balance matrices coincide with the covariates used by the corresponding
analysis specification: cubic or saturated age terms in Vietnam, the Card or
Kitagawa conditioning set in the schooling application, and the common six
covariates in the Child application. Only algebraic or numerically verified
outcome-weight vectors enter these plots. Because the adjusted weights can be
negative, a reported binary weighted mean is an algebraic signed moment and
need not be a probability or lie in the unit interval. The plots therefore
describe balance in the included observed covariate means; they do not test
instrument validity, balance unobserved covariates, or identify the source of
any remaining difference.

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

All six application-specific XGBoost tuning notebooks—Wald-AIPW and PLR-IV for
Vietnam, Card, and Child—default to loading their verified exports with
`RERUN_TUNING = FALSE`. Setting this switch to `TRUE` deliberately repeats the
expensive fold-specific nested-tuning procedure and replaces the corresponding
tuning export after the run has completed.

Complete-pipeline transformation checks rerun the relevant nuisance estimation,
cross-fitting, and tuning rules after the outcome has been shifted. They are
distinct from frozen-weight algebraic checks, which hold a fitted weight vector
fixed and evaluate the implied change directly.

The three headline translation notebooks classify a complete-pipeline rerun as
unchanged when the absolute difference between the separately fitted original-
and shifted-outcome estimates does not exceed $10^{-4}$. This is an operational
reporting rule for the empirical rerun. It is separate from both R's
`all.equal()` outcome-weight reconstruction check and the $10^{-4}$ display
epsilon used when presenting small individual weights.

The supplementary Child four-unit notebook studies maternal log income in
dollars, cents, thousands of dollars, and hundred-thousands of dollars. It uses
the stricter absolute threshold $10^{-8}$ by design to inspect whether the
translation result continues to hold numerically under much larger positive
and negative log-unit changes. This stricter diagnostic does not replace the
$10^{-4}$ reporting convention used by the headline Child translation table.

The shared `translation_rerun_row()` helper now defaults to the headline
$10^{-4}$ reporting tolerance, and every headline translation notebook also
passes that tolerance explicitly. The supplementary Child four-unit notebook
does not rely on this helper default: it applies its deliberately stricter
$10^{-8}$ diagnostic rule explicitly. Callers should continue to pass their
intended tolerance so that the applicable reporting convention remains visible
at the point where results are classified.

### Translation XGBoost runtime switch

Each main translation notebook defaults to `RUN_XGB_RERUN = FALSE`. In this
state, the notebook loads a small verified `*_translation_xgb_rerun_export.rds`
file, checks its schema, required estimator/specification rows, seed, fold
counts, tuning budget, and `tune_on_folds` setting, and then recalculates the
$10^{-4}$ classification from the stored unrounded rerun differences. The
XGBoost translation table is therefore complete without access to a local
knitr cache. Neither the XGBoost translation export nor the combined complete
translation export is overwritten in this default state.

Setting `RUN_XGB_RERUN = TRUE` is a deliberate rebuild. It requires the two
matching original-outcome XGBoost tuning exports, repeats shifted-outcome
nested tuning in a non-cached chunk, validates the rebuilt rows, and replaces
the application-specific XGBoost and complete translation exports. The final
save is refused when an expected result family is missing. This switch affects
the shifted-outcome XGBoost branch only; the other translation sections retain
their existing individual runtime and cache settings.

For the supplementary Child four-unit diagnostic, each learner receives an
outcome constructed directly as `Y_dollars + shift`. Expressions based on
rescaling the original income before taking logs are evaluated only as a
numerical unit-coding check. This avoids feeding adaptive learners outcomes
that are mathematically equivalent but differ slightly because of floating-
point evaluation order.

## Reproducibility

The repository-level `README.md` records the validated R and package versions,
required data files, notebook execution order, stored tuning exports, and cache
policy. `REPRODUCING.md` provides the corresponding step-by-step instructions
for a fresh checkout without local caches. For the final replication
repository, retain the canonical `functions_all.R`, the application notebooks,
the verified `.rds` exports, and a repository-level `renv.lock`. Do not
distribute duplicate helper files inside the application folders.

The stored artifacts have distinct roles. The six `*_xgb*_tuning_export.rds`
files supply the original-outcome XGBoost fits and extracted weights. The three
`*_translation_xgb_rerun_export.rds` files supply verified shifted-outcome
XGBoost rows during an ordinary translation render. The three
`*_translation_rerun_results.rds` files are the combined translation
references. Knitr caches are disposable runtime accelerators, not replication
inputs. The older Card and Child `*_translation_xgb_method_a_export.rds` files
are not read by the current pipeline and are superseded by the corresponding
`*_translation_xgb_rerun_export.rds` files.

Seeds, fold counts, and stored exports control the current analysis, but exact
cross-machine equality of forest and XGBoost reruns is not yet guaranteed.
Package versions, thread scheduling, and cached objects can affect fitted
learners. A final clean replication should therefore pin the R environment,
use explicit single-thread settings where exact equality is required, and run
the validation once without relying on pre-existing knitr caches.

On 4 August 2026, the three headline translation notebooks were rendered after
moving aside their relevant caches. The resulting 28 Vietnam, 56 Card, and 28
Child rows matched the stored complete-result references exactly. The Child
four-unit notebook was then rerun from a new cache after constructing every
income coding as an exact additive shift. All 28 dollar/cent comparisons with
Method A passed its $10^{-8}$ audit, with a maximum absolute difference of
$4.76\times10^{-13}$. This resolves the previously observed cross-notebook
honest-forest and Ranger discrepancy on the validated machine.

The application-specific smoother-condition audits cover the same eligible
DML learner families used in the outcome-weight diagnostics. The Card audit
covers its four conditioning-set/treatment cells, the Vietnam audit covers its
cubic and saturated age specifications, and the Child audit covers the work-
status and maternal-log-income samples. XGBoost vectors that fail the estimator
reproduction gate are not given a smoother-condition interpretation.

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
