# Vietnam notebooks — script map

 This
file is a local map of the notebooks and the order in which to run them. For the
R environment, runtime switches, data checksums, and rebuild instructions, see
the repository-root [`README.md`](../README.md) and
[`REPRODUCING.md`](../REPRODUCING.md).

Shared estimator and diagnostic functions come from the repository-root
`functions_all.R`, loaded with `source(file.path("..", "functions_all.R"))`.

## Notebooks

| Notebook | Purpose | Thesis output | Reads / writes |
|---|---|---|---|
| `vietnam_xgb_tuning.Rmd` | Nested XGBoost tuning for the Wald-AIPW nuisances | DML learner comparison | reads `sipp.dta`; reads/writes `vietnam_xgb_tuning_export.rds` |
| `vietnam_xgb_pliv_tuning.Rmd` | Nested XGBoost tuning for the PLR-IV nuisances | DML learner comparison | reads `sipp.dta`; reads/writes `vietnam_xgb_pliv_tuning_export.rds` |
| `vietnam_presentation_4.Rmd` | Headline kappa-vs-DML point estimates and outcome-weight normalization, concentration, and balance | Section 5 point-estimate and normalization tables | reads `sipp.dta` and the two tuning exports |
| `vietnam_translation_invariance_method_a.Rmd` | Complete-pipeline translation reruns (kappa, honest forest, OLS/ranger, XGBoost) | Translation-invariance table | reads `sipp.dta`, `vietnam_translation_xgb_rerun_export.rds`; writes `vietnam_translation_rerun_results.rds` on a deliberate rebuild only |
| `vietnam_descriptives_1.Rmd` | Sample construction, first stage, propensity overlap | Appendix design diagnostics and propensity figure | reads `sipp.dta` |
| `vietnam_smoother_conditions_ols_ranger.Rmd` | Condition 3 / 5 smoother audit for the OLS and ranger learners | Appendix | reads `sipp.dta` |
| `vietnam_forest_translation_smoother_audit.Rmd` | Localizes the saturated honest-forest translation departure to the outcome smoother | Appendix smoother-change table | reads `sipp.dta` |

## Run order

1. `vietnam_xgb_tuning.Rmd` with `RERUN_TUNING = FALSE` (loads the verified export)
2. `vietnam_xgb_pliv_tuning.Rmd`  with `RERUN_TUNING = FALSE`
3. `vietnam_presentation_4.Rmd`
4. `vietnam_translation_invariance_method_a.Rmd`  with `RUN_XGB_RERUN = FALSE`
5. `vietnam_descriptives_1.Rmd`
6. `vietnam_smoother_conditions_ols_ranger.Rmd` —optional appendix audit
7. `vietnam_forest_translation_smoother_audit.Rmd`  optional internal audit

The two tuning notebooks come first only so their verified `.rds` exports are
present for the presentation and translation notebooks. Under the default
switches they load those exports rather than repeating the expensive tuning.
