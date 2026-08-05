# Outcome-Weight Diagnostics for Instrumental-Variables Estimators

This repository contains the empirical code for a thesis comparing classical
kappa-weighting estimators with double/debiased machine-learning estimators.
The analysis studies point estimates and their implied outcome weights across
three instrumental-variables applications:

1. Vietnam-era draft eligibility, military service, and wages;
2. proximity to college, schooling, and wages (Card);
3. sibling-sex composition, fertility, and maternal labor and income outcomes
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

A beginner-facing checklist for a fresh checkout, required data and exports,
default switches, ordinary rendering, and deliberate rebuilds is provided in
[`REPRODUCING.md`](REPRODUCING.md).

## Repository structure

```text
.
├── functions_all.R       # canonical shared estimator/diagnostic library
├── README.md             # project overview
├── IMPLEMENTATION.md     # estimator and diagnostic details
├── REPRODUCING.md        # clean-checkout and rerun instructions
├── vietnam_final/        # Vietnam application
├── card_final/           # Card application
└── child/                # Angrist–Evans application
```

The working thesis directory may additionally contain the LaTeX source,
bibliography, figures, and compiled PDF. They are not required in the empirical
replication repository unless the thesis itself is also being distributed.

`functions_all.R` is the only active shared function file. Application
notebooks load it with:

```r
source(file.path("..", "functions_all.R"))
```

Do not place additional copies inside the application directories. Keeping one
canonical helper prevents estimator definitions and diagnostics from drifting
across applications.

## Software environment

The analysis was last validated on 4 August 2026 with R 4.5.3 and the following
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

The current working folder does not yet contain a complete `renv` environment.
Until one is created, the table above documents the validated environment but
does not automatically restore it. After the final repository has been
assembled, initialize `renv` once from that repository's root:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

Commit `renv.lock`, `.Rprofile`, `renv/activate.R`, and
`renv/settings.json` if it is created. Do not commit `renv/library/`.
These files form one repository-level environment; the three application
folders do not need separate lockfiles. Never distribute an `.Rprofile` that
refers to `renv/activate.R` unless that activation file is present as well.
The complete author and professor workflows are given in `REPRODUCING.md`.

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

All six XGBoost tuning notebooks—Wald-AIPW and PLR-IV for Vietnam, Card, and
Child—default to:

```r
RERUN_TUNING <- FALSE
```

In this mode they read the verified `.rds` exports rather than repeating
expensive nested tuning. The main application notebooks also import these
exports. Keep the exports when an ordinary, reasonably fast render is desired.

The three main translation notebooks default to:

```r
RUN_XGB_RERUN <- FALSE
```

In this mode they load and validate the verified shifted-outcome XGBoost
translation exports. No XGBoost cache, tuning, or weight extraction is needed,
and the official translation result files are not overwritten. Setting the
switch to `TRUE` deliberately repeats shifted-outcome nested tuning with the
documented protocol. The XGBoost rerun chunk is then evaluated without a cache,
validated, and saved together with the complete translation result.

## Recommended execution order

### Vietnam

1. `vietnam_xgb_tuning.Rmd`
2. `vietnam_xgb_pliv_tuning.Rmd`
3. `vietnam_presentation_4.Rmd`
4. `vietnam_translation_invariance_method_a.Rmd`
5. `vietnam_descriptives_1.Rmd`
6. `vietnam_smoother_conditions_ols_ranger.Rmd`, when the additional smoother
   audit is required

The ordinary tuning-notebook render loads:

```text
vietnam_xgb_tuning_export.rds
vietnam_xgb_pliv_tuning_export.rds
```

The translation notebook combines newly fitted non-XGBoost rows with its
verified XGBoost rows in:

```text
vietnam_translation_rerun_results.rds
vietnam_translation_xgb_rerun_export.rds
```

With the default `RUN_XGB_RERUN = FALSE`, both files are treated as verified
references and are not replaced. With `TRUE`, the shifted-outcome XGBoost
branch is rebuilt and the validated files are replaced.

### Card

1. `card_xgb_tuning.Rmd`
2. `card_xgb_pliv_tuning.Rmd`
3. `card_test.Rmd`
4. `card_translation_invariance_method_a_2.Rmd`
5. `card_xgb_sensitivity.Rmd`
6. `card_descriptives_1.Rmd`
7. `card_smoother_conditions.Rmd`, when the additional smoother audit is
   required
8. Run `card_final/export_thesis_figures.R` from the repository root when the
   thesis figures must be refreshed

The Card tuning notebooks read or rebuild
`card_xgb_tuning_export.rds` and `card_xgb_pliv_tuning_export.rds`. The
translation notebook writes `card_translation_rerun_results.rds` and
`card_translation_xgb_rerun_export.rds`; the additional smoother audit writes
`card_smoother_condition_results.rds`.

### Child

1. `ae98_make_balanced_subsample_submission.Rmd`
2. `child_xgb_tuning.Rmd`
3. `child_xgb_pliv_tuning.Rmd`
4. `child_test.Rmd`
5. `child_translation_invariance_method_a.Rmd`
6. `ae98_descriptives.Rmd`
7. `child_smoother_conditions.Rmd`, when the smoother-condition audit is
   required
8. `child_income_translation_four_units.Rmd`, only for the supplementary
   strict unit-change diagnostic

The sample-construction notebook writes `ae98_sub.rds`. The two tuning
notebooks read or rebuild `ae98_xgb_tuning_export.rds` and
`ae98_xgb_pliv_tuning_export.rds`. The headline translation notebook writes
`ae98_translation_rerun_results.rds` and
`ae98_translation_xgb_rerun_export.rds`. The supplementary four-unit notebook
writes separate `ae98_income_four_unit_*.rds` files and does not overwrite the
headline translation or tuning exports. The Child smoother-condition audit
writes `ae98_smoother_condition_results.rds`.

## Repeating the XGBoost tuning

For any application, deliberately repeat a tuning exercise by changing the
corresponding Wald-AIPW or PLR-IV tuning notebook to:

```r
RERUN_TUNING <- TRUE
```

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

## Numerical comparison conventions

The code intentionally uses different numerical rules for different tasks.
They should not be replaced by one common threshold.

- Outcome-weight reconstruction is checked by `check_weight_identity()`, which
  calls R's `all.equal()` without overriding its default scale-aware tolerance.
- The three headline complete-pipeline translation notebooks classify an
  estimate as unchanged when the absolute rerun difference is at most
  $10^{-4}$.
- `child_income_translation_four_units.Rmd` is a supplementary strict
  diagnostic. It compares maternal log income measured in dollars, cents,
  thousands, and hundred-thousands and deliberately uses $10^{-8}$.
- Propensity-score clipping at $10^{-8}$ and tighter checks used for sample
  metadata, optimizers, and matrix calculations are numerical safeguards, not
  translation-invariance reporting rules.

The shared `translation_rerun_row()` function defaults to $10^{-4}$, and all
headline translation notebooks also pass $10^{-4}$ explicitly. The strict
Child four-unit diagnostic applies its separate $10^{-8}$ rule explicitly.

## Stored results, caches, and generated files

The `.rds` exports are small, intentional pipeline artifacts. Downstream
notebooks read them, and they avoid repeating the most expensive fits.

Two older files,
`card_final/card_translation_xgb_method_a_export.rds` and
`child/ae98_translation_xgb_method_a_export.rds`, are not read by the current
pipeline. They are superseded by the corresponding
`*_translation_xgb_rerun_export.rds` files and should be omitted from the final
replication repository.

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

For the final replication audit, first preserve the current files and exports,
then run the selected validation in a fresh R session without using the
affected notebook cache. Cache directories should not be deleted before the
current state has been versioned or backed up. In particular, changing
`functions_all.R` does not by itself guarantee that every previously cached
model fit will be recomputed.

## Current reproducibility limitations

The stored exports make the present results inspectable and avoid accidental
retuning, but the repository is not yet a fully restorable clean-machine
package. It currently has no `renv.lock`, the fitted forest and XGBoost
learners do not impose a common explicit single-thread setting, and result
exports do not embed checksums of the raw data, notebook, and shared helper.
The validated raw-data checksums are now recorded in `REPRODUCING.md`. The final
external test must use a fresh clone on a second machine after the environment
lockfile has been added.

The earlier supplementary Child discrepancy has been resolved. The four-unit
notebook now constructs every unit as an exact additive shift of the same
dollar-log vector. In a cache-free run on 4 August 2026, its 28 dollar/cent
comparisons all reproduced the corresponding Method A estimates under the
strict $10^{-8}$ audit rule. The maximum absolute cross-notebook difference was
$4.76\times10^{-13}$. The cache-free headline Child translation notebook also
reproduced all 28 stored result rows exactly. No additional thread control was
needed for these results on the validated machine.

## Validation status

At the 4 August 2026 static audit:

- all 21 R Markdown notebooks, including the new Child smoother-condition
  audit, parsed successfully;
- `functions_all.R` and both standalone Card scripts parsed successfully;
- every notebook that uses shared estimator functions referenced the root
  helper correctly;
- the six kappa outcome-weight identities held on the three empirical datasets
  to numerical precision;
- CBPS converged for all audited specifications;
- the Vietnam main presentation completed all 74 chunks;
- cache-free headline translation renders reproduced all 28 Vietnam, 56 Card,
  and 28 Child stored rows, with no numerical or classification discrepancies;
- the cache-free Child four-unit diagnostic reproduced all 28 corresponding
  Method A dollar/cent rows and generated all 56 four-unit rows;
- the Child smoother-condition notebook produced its 18 Condition 3 and 18
  Condition 5 assessments, and every included outcome-weight vector reproduced
  its fitted estimate;
- the refreshed translation estimates agreed with the thesis at the reported
  precision.

The XGBoost reconstructed outcome weights do not exactly reproduce their fitted
estimates under the identity check. The notebooks retain those fitted estimates
for learner comparisons but exclude failed reconstructed vectors from
weight-based and balance diagnostics.

## Suggested GitHub workflow

Use one R project and one `renv` environment at the repository root, and
preserve the directory structure shown above. After the final lockfile and
authorized data have been supplied, a collaborator can start a fresh R session
in that root and run:

```r
install.packages("renv")
renv::restore(prompt = FALSE)
renv::status()
rmarkdown::render("vietnam_final/vietnam_presentation_4.Rmd")
```

The committed root `.Rprofile` loads `renv/activate.R` automatically when R is
started in the repository. The collaborator does not source the lockfile or
activation script manually. `renv::restore()` installs the versions recorded
in `renv.lock`; it does not run the empirical notebooks.

For the clearest RStudio workflow, include one `.Rproj` file in the repository
root and omit the three old application-level `repro.Rproj` files. Opening an
application-level project would make that subfolder the project root and would
not activate the repository-level environment in the intended way.

Do not upload local knitr caches or `renv/library/`. Do upload the source
notebooks, `functions_all.R`, documentation, the lockfile, and any datasets or
verified result exports that redistribution rules permit.
