# Current thesis structure

Status: synchronized with `draft_of_13.tex` on 5 August 2026. This is an
internal writing outline, not a replication instruction. Omit it from the final
empirical-code repository unless a thesis outline is intentionally included.
The authoritative manuscript is always the LaTeX source.

## Research question

How do observation-level outcome weights characterize the finite-sample
behavior of kappa weighting estimators and data-adaptive IV estimators?

The empirical analysis addresses three subsidiary questions:

1. Do the fitted estimators satisfy the relevant normalization and
   outcome-invariance properties?
2. Do estimators with similar point estimates nevertheless differ in effective
   support, weight concentration, sign cancellation, or covariate balance?
3. How sensitive are these properties to the instrument-propensity
   specification, nuisance learner and tuning choices, and empirical IV design?

## Main text

### 1. Introduction

- Motivation: economically irrelevant outcome recodings should not alter the
  estimated effect.
- Related literature: conditional LATE identification, kappa weighting,
  normalization, DML, and outcome-weight analysis.
- Research question and empirical approach.
- Contribution: compare five kappa estimator forms, reported as six
  specifications because the normalized Uysal estimator is fitted with MLE and
  CBPS instrument propensities, with PLR-IV and Wald-AIPW implementations.

### 2. Econometric Framework

- Potential outcomes, compliance types, and the conditional LATE assumptions.
- Abadie's kappa representation of complier moments.
- The distinction between kappa identification weights and realized outcome
  weights.

### 3. Kappa and DML-IV estimators

- Five kappa-based LATE estimator forms and the normalized versus unnormalized
  distinction.
- Maximum-likelihood and covariate-balancing instrument propensity scores.
- PLR-IV as a partially linear IV estimator; it should not automatically be
  described as identifying the same LATE without the required structural
  assumptions.
- Wald-AIPW as the DML LATE score.
- Linear/logistic, honest-forest, Ranger, and XGBoost nuisance learners.
- Five outer folds, fold-specific XGBoost tuning with three inner folds and 15
  random-search evaluations.

### 4. Outcome Weights, Invariance, and Empirical Diagnostics

- Pseudo-IV representation and the outcome-weight identity.
- Scale-normalized and fully normalized outcome-weight classes.
- Closed-form outcome weights of the normalized Uysal estimator.
- Translation and logarithmic scale invariance.
- Complete-pipeline reruns after outcome transformations, distinguished from a
  frozen-weight algebraic check.
- Normalization, effective support, concentration, sign composition,
  cancellation, and covariate-balance diagnostics.

### 5. Vietnam Draft Lottery

- Point estimates under alternative age specifications.
- Complete-pipeline translation invariance.
- Outcome-weight normalization, support, concentration, and balance.
- Additional smoother-condition results remain in the appendix.

### 6. College Proximity and Returns to Schooling

- Two treatment margins and two conditioning sets.
- Point estimates and learner sensitivity.
- Complete-pipeline translation invariance.
- Outcome-weight normalization, concentration, and balance.
- XGBoost sensitivity and smoother-condition details remain in the appendix.

### 7. Fertility and Mothers' Labor-Market Outcomes

- Fixed labor-status and positive-income samples of 3,000 observations each.
- Point estimates and nuisance-learner sensitivity.
- Complete-pipeline translation invariance for labor status and maternal log
  income.
- Supplementary four-unit log-income diagnostic.
- Outcome-weight normalization, concentration, and balance.

### 8. Discussion

- Cross-application differences in assignment, first-stage strength,
  conditioning complexity, and sample construction.
- Distinguish realized outcome-weight normalization from complete-pipeline
  invariance.
- Interpret effective support, concentration, sign cancellation, and balance as
  separate diagnostics, not as a single estimator ranking.
- Discuss learner-specific results and the limitation that the extracted
  XGBoost weights fail the estimate-reproduction gate.
- Limitations: diagnostic rather than identifying content, implementation
  specificity, fixed Child subsamples, positive-income selection, and the
  XGBoost weight-extraction limitation.

### 9. Conclusion

- Outcome weights make the realized empirical comparison explicit.
- Complete-pipeline reruns reveal whether the implemented estimator remains
  stable under economically irrelevant outcome transformations.
- Similar coefficients do not imply similar finite-sample weighting behavior.

## Appendices

- Algebraic normalization of the kappa estimators.
- PLR-IV as an aggregation of conditional Wald estimands.
- Outcome-weight classification and kappa derivations.
- Full Vietnam, Card, and Child tables and diagnostics.
- Instrument-propensity, first-stage, overlap, tuning, translation, smoother,
  concentration, and covariate-balance details.

## Terminology controls

- Say **five kappa estimator forms** and **six reported kappa specifications**.
- Say **normalized Uysal estimator** for `tau_u`; `tau_a,10` is the other
  normalized form and should not be called the Uysal estimator.
- Use **normalized-form** and **unnormalized-form** when grouping kappa
  estimators.
- Use **scale-normalized** and **fully normalized** only for Knaus's realized
  outcome-weight classes.
- Keep PLR-IV and Wald-AIPW as separate estimator families.
- Do not interpret XGBoost normalization, support, concentration, sign, or
  balance when its reconstructed weights fail to reproduce the fitted estimate.
