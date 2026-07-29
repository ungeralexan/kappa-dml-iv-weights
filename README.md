# Outcome-Weight Diagnostics for Instrumental-Variables Estimators

This repository contains the empirical code for a thesis comparing classical
kappa-weighting estimators with double/debiased machine-learning estimators.
The analysis studies point estimates and their implied outcome weights across
three instrumental-variables applications:

1. Vietnam-era draft eligibility, military service, and wages;
2. proximity to college, schooling, and wages (Card);
3. sibling-sex composition, fertility, and maternal labor outcomes
   (Angrist–Evans).

The code compares 2SLS, normalized and unnormalized kappa estimators, PLR-IV,
and Wald-AIPW under linear/logistic, honest-forest, Ranger, and XGBoost
nuisance learners. Diagnostics cover numerical weight identities,
normalization, effective sample size, sign cancellation, covariate balance,
learner sensitivity, and translation invariance.

Detailed notes on the custom kappa implementation, propensity-score fitting,
analytical standard errors, outcome-weight reconstruction, and its relationship
to available R and Stata software are provided in
[`IMPLEMENTATION.md`](IMPLEMENTATION.md).

## Repository structure

```text
.
├── functions_all.R       # canonical shared estimator/diagnostic library
├── vietnam_final/        # Vietnam application
├── card_final/           # Card application
├── child/                # Angrist–Evans application
├── draft_of_13.tex       # thesis source
├── references.bib
├── pics/                 # figures used by the thesis
└── output/pdf/           # compiled thesis output
```

`functions_all.R` is the only active shared function file. Application
notebooks load it with:

```r
source(file.path("..", "functions_all.R"))
```

Do not place additional copies inside the application directories. Keeping one
canonical helper prevents estimator definitions and diagnostics from drifting
across applications.

## Software environment

The analysis was last validated on 29 July 2026 with R 4.5.3 and the following
package versions:

| Package | Version |
|---|---:|
| AER | 1.2.16 |
| cobalt | 4.6.2 |
| dplyr | 1.2.1 |
| DoubleML | 1.0.2 |
| ggplot2 | 4.0.2 |
| grf | 2.6.1 |
| gridExtra | 2.3 |
| haven | 2.5.5 |
| kableExtra | 1.4.0 |
| knitr | 1.51 |
| lmtest | 0.9.40 |
| mlr3 | 1.6.0 |
| mlr3learners | 0.14.0 |
| mlr3tuning | 1.6.0 |
| mlr3verse | 0.3.1 |
| OutcomeWeights | 0.2.0.9000 |
| ranger | 0.18.0 |
| rmarkdown | 2.31 |
| sandwich | 3.1.1 |
| scales | 1.4.0 |
| tidyr | 1.3.2 |
| viridis | 0.6.5 |
| xgboost | 3.2.1.1 |

The installed `OutcomeWeights` version came from
`MCKnaus/OutcomeWeights` at commit
`0c94f940b04c14d0247b46842af37752e306b79e`.

The project does not yet contain an `renv.lock`. Until one is created, the
table above documents the validated environment but does not automatically
restore it. To initialize `renv` from the repository root:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

Commit `renv.lock`, `.Rprofile`, `renv/activate.R`, and
`renv/settings.json` if it is created. Do not commit `renv/library/`.

## Data inputs

The application notebooks currently expect:

```text
vietnam_final/sipp.dta
card_final/card.dta
child/ae98.dta
```

The Child analysis first creates `child/ae98_sub.rds`, containing the fixed
labor and positive-income analysis samples. Before publishing the repository,
confirm that the three source datasets may legally be redistributed. If not,
replace them with download and preparation instructions.

## Ordinary rendering

Run a notebook from its own application directory, or use
`rmarkdown::render()` from the repository root. For example:

```r
rmarkdown::render("vietnam_final/vietnam_presentation_4.Rmd")
```

The Vietnam and Card XGBoost tuning notebooks default to:

```r
RERUN_TUNING <- FALSE
```

In this mode they read the verified `.rds` exports rather than repeating
expensive nested tuning. The main application notebooks also import these
exports. Keep the exports when an ordinary, reasonably fast render is desired.

The two Child tuning notebooks currently do not expose a
`RERUN_TUNING` switch and rerun their tuning exercises when rendered. Their
verified exports are nevertheless read by `child_test.Rmd` and the downstream
translation analysis. Avoid rendering the Child tuning notebooks unless a
deliberate retuning run is intended.

## Recommended execution order

### Vietnam

1. `vietnam_xgb_tuning.Rmd`
2. `vietnam_xgb_pliv_tuning.Rmd`
3. `vietnam_presentation_4.Rmd`
4. `vietnam_translation_invariance_method_a.Rmd`
5. `vietnam_descriptives_1.Rmd`
6. `vietnam_smoother_matrix_diagnostics.Rmd`, when the additional smoother
   audit is required

The ordinary tuning-notebook render loads:

```text
vietnam_xgb_tuning_export.rds
vietnam_xgb_pliv_tuning_export.rds
```

The translation notebook performs a complete shifted-outcome rerun and writes:

```text
vietnam_translation_rerun_results.rds
vietnam_translation_xgb_rerun_export.rds
```

### Card

1. `card_xgb_tuning.Rmd`
2. `card_xgb_pliv_tuning.Rmd`
3. `card_test.Rmd`
4. `card_translation_invariance_method_a_2.Rmd`
5. `card_xgb_sensitivity.Rmd`
6. `card_descriptives_1.Rmd`
7. Run `card_final/export_thesis_figures.R` from the repository root

### Child

1. `ae98_make_balanced_subsample_submission.Rmd`
2. `child_xgb_tuning.Rmd`
3. `child_xgb_pliv_tuning.Rmd`
4. `child_test.Rmd`
5. `child_translation_invariance_method_a.Rmd`
6. `ae98_descriptives.Rmd`

## Repeating the XGBoost tuning

For Vietnam and Card, deliberately repeat a tuning exercise by changing the
corresponding notebook to:

```r
RERUN_TUNING <- TRUE
```

The Child tuning notebooks currently retune whenever they are rendered.

The tuning design uses seed 42, five outer cross-fitting folds, three inner
folds, fifteen random-search evaluations, and fold-specific tuning. These
settings make reruns controlled, but exact cross-machine equality is not yet
guaranteed. Package versions, R versions, operating systems, hardware, and
multithreaded XGBoost calculations can produce small differences and may
occasionally select different hyperparameters.

For the most defensible replication package:

1. create and commit `renv.lock`;
2. retain the verified tuning exports;
3. record the R version and platform;
4. use explicit single-thread XGBoost settings if bit-level reproducibility is
   required;
5. distinguish ordinary rendering from deliberate full retuning.

## Stored results, caches, and generated files

The `.rds` exports are small, intentional pipeline artifacts. Downstream
notebooks read them, and they avoid repeating the most expensive fits.

Knitr directories such as `*_cache/`, `cache/`, and `*_files/` are generated
locally. They are ignored by Git. Deleting a cache does not delete the source
analysis, but it forces the affected chunks to be recomputed. Figure folders
may also be required by the LaTeX thesis until the relevant notebook has been
rendered again.

Rendered `.html` files, R session files (`.RData`, `.Rhistory`), `.DS_Store`,
RStudio project state, temporary audit files, and LaTeX auxiliary files are not
source inputs and are excluded by `.gitignore`.

## Cache-related errors

After changing `functions_all.R`, restart R before rendering:

- RStudio: **Session → Restart R**
- macOS shortcut: **Cmd + Shift + F10**

If a notebook then reports a table-schema or missing-object error, its knitr
cache may contain objects produced by the previous helper version. Remove or
rename only that notebook's cache directory and render it again. Do not remove
the tuning `.rds` exports unless the intention is to repeat tuning.

## Validation status

At the last audit:

- all 18 R Markdown notebooks parsed successfully;
- all standalone R files parsed successfully;
- every notebook that uses shared estimator functions referenced the root
  helper correctly;
- the six kappa outcome-weight identities held on the three empirical datasets
  to numerical precision;
- CBPS converged for all audited specifications;
- the Vietnam main presentation completed all 74 chunks;
- the refreshed Vietnam estimates agreed with the thesis at the reported
  precision.

The XGBoost reconstructed outcome weights do not exactly reproduce their fitted
estimates under the identity check. The notebooks retain those fitted estimates
for learner comparisons but exclude failed reconstructed vectors from
weight-based and balance diagnostics.

## Suggested GitHub workflow

Use the repository root as the R project root and preserve the directory
structure shown above. A collaborator can then:

```r
install.packages("renv")
renv::restore()  # after renv.lock has been committed
rmarkdown::render("vietnam_final/vietnam_presentation_4.Rmd")
```

Do not upload local knitr caches or `renv/library/`. Do upload the source
notebooks, `functions_all.R`, documentation, the lockfile, and any datasets or
verified result exports that redistribution rules permit.
