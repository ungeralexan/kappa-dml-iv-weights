# Reproducing the thesis analysis

This guide is written for a reader who has downloaded or cloned the project
without any of the author's local knitr cache directories. Those caches are not
required. The small verified `.rds` exports are intentional pipeline inputs and
must remain in the project.

## 1. Required software

The analysis was validated with R 4.5.3. RStudio is convenient but not
required. The exact package versions currently used are listed in `README.md`.
The project does not yet contain an `renv.lock`, so `renv::restore()` is not yet
an available installation route. Until the lockfile is added, install the
listed package versions manually, including the development version of
`OutcomeWeights` from commit
`0c94f940b04c14d0247b46842af37752e306b79e`.

The lockfile will be added as a separate reproducibility step after the current
notebook workflow is stable. When present, there will be one `renv.lock` at the
repository root, not one in every application folder.

## 2. Required data layout

Place the raw datasets at exactly these relative paths:

```text
vietnam_final/sipp.dta
card_final/card.dta
child/ae98.dta
```

The Child workflow also uses the fixed subsample export:

```text
child/ae98_sub.rds
```

If `ae98_sub.rds` is unavailable, render
`child/ae98_make_balanced_subsample_submission.Rmd` first. It uses the fixed
documented sampling seed and refuses to replace an existing fixed sample when
the regenerated row identifiers differ.

Do not change the notebooks to point to a personal Dropbox or home-directory
path. Download the datasets from the authorized location and place them in the
folders above.

## 3. Required verified exports

Ordinary rendering does not repeat expensive XGBoost tuning or weight
extraction. Keep these tuning exports:

```text
vietnam_final/vietnam_xgb_tuning_export.rds
vietnam_final/vietnam_xgb_pliv_tuning_export.rds
card_final/card_xgb_tuning_export.rds
card_final/card_xgb_pliv_tuning_export.rds
child/ae98_xgb_tuning_export.rds
child/ae98_xgb_pliv_tuning_export.rds
```

The main translation notebooks also load these verified complete-rerun XGBoost
exports by default:

```text
vietnam_final/vietnam_translation_xgb_rerun_export.rds
card_final/card_translation_xgb_rerun_export.rds
child/ae98_translation_xgb_rerun_export.rds
```

The complete combined translation references are:

```text
vietnam_final/vietnam_translation_rerun_results.rds
card_final/card_translation_rerun_results.rds
child/ae98_translation_rerun_results.rds
```

These files are small and should be distributed with the code. Knitr cache
folders such as `*_cache/` and `child/cache/` are large local speed-ups and
should not be distributed.

## 4. Default switches

The six original-outcome XGBoost tuning notebooks use:

```r
RERUN_TUNING <- FALSE
```

With this default, the notebook reads its verified tuning export, including
the previously extracted outcome weights. Set this switch to `TRUE` only to
repeat the expensive original-outcome nested tuning and weight extraction.

The three main translation notebooks use:

```r
RUN_XGB_RERUN <- FALSE
```

With this default, the notebook loads and validates its verified shifted-
outcome XGBoost translation export. It does not need an XGBoost cache and does
not replace either the XGBoost translation export or the complete combined
translation result.

Setting `RUN_XGB_RERUN <- TRUE` deliberately retunes the shifted-outcome
XGBoost models. The XGBoost chunk has caching disabled, so `TRUE` means a fresh
rerun. After schema and setting checks pass, the notebook replaces both its
XGBoost translation export and its complete combined translation result.

## 5. Ordinary translation render without caches

From the repository root, start a fresh R session and run:

```r
rmarkdown::render("vietnam_final/vietnam_translation_invariance_method_a.Rmd")
rmarkdown::render("card_final/card_translation_invariance_method_a_2.Rmd")
rmarkdown::render("child/child_translation_invariance_method_a.Rmd")
```

Under the default switches, kappa and the non-XGBoost learner sections are
evaluated according to their existing chunk-cache settings. On a new computer
their local caches do not exist, so those sections are fitted and new caches
are created automatically. The expensive XGBoost translation rows are read
from the verified exports. Official result exports remain unchanged.

The first uncached render can therefore take longer than later renders even
though XGBoost is not retuned.

## 6. Deliberate XGBoost translation rebuild

Before rebuilding, preserve the committed reference files. Then change only
the relevant translation notebook to:

```r
RUN_XGB_RERUN <- TRUE
```

Render that notebook. A rebuild requires the corresponding two original-
outcome tuning exports. It uses seed 42, five outer folds, three inner folds,
15 random-search evaluations, and `tune_on_folds = TRUE`. The newly generated
XGBoost rows are checked before the official translation exports are replaced.

After inspecting the new results, return the switch to:

```r
RUN_XGB_RERUN <- FALSE
```

## 7. Child four-unit diagnostic

`child/child_income_translation_four_units.Rmd` is a separate strict diagnostic.
It constructs cents, thousands, and hundred-thousands as exact additive shifts
of the same dollar-log vector and uses an absolute threshold of $10^{-8}$. Its
learner switches currently request substantive reruns, including XGBoost, so
rendering it can be expensive. It writes separate `ae98_income_four_unit_*.rds`
files and does not replace the headline translation or tuning exports.

The cache-free validation on 4 August 2026 generated 56 four-unit rows. All 28
dollar/cent comparisons against the stored Method A results passed, with a
maximum absolute difference of $4.76\times10^{-13}$. The current output files
are:

```text
child/ae98_income_four_unit_translation_results.rds
child/ae98_income_four_unit_xgb_export.rds
```

## 8. Recommended application order

Vietnam:

1. Render the two XGBoost tuning notebooks with `RERUN_TUNING = FALSE` to check
   that their verified exports load.
2. Render `vietnam_presentation_4.Rmd`.
3. Render `vietnam_translation_invariance_method_a.Rmd`.
4. Render the descriptives and optional smoother audit.

Card:

1. Render the two XGBoost tuning notebooks with `RERUN_TUNING = FALSE`.
2. Render `card_test.Rmd`.
3. Render `card_translation_invariance_method_a_2.Rmd`.
4. Render the sensitivity, descriptives, and smoother notebooks.
5. Run `card_final/export_thesis_figures.R` from the repository root when the
   thesis figures must be refreshed.

Child:

1. Verify or create `ae98_sub.rds`.
2. Render the two XGBoost tuning notebooks with `RERUN_TUNING = FALSE`.
3. Render `child_test.Rmd`.
4. Render `child_translation_invariance_method_a.Rmd`.
5. Render descriptives.
6. Render `child_smoother_conditions.Rmd` when the smoother-condition audit is
   required. It writes `ae98_smoother_condition_results.rds`.
7. Run the strict four-unit notebook only when that additional diagnostic is
   deliberately requested.

## 9. What reproduction means

The replication target is agreement of samples, reported estimates, standard
errors, and translation classifications at the precision used in the thesis.
Binary-identical cache files or `.rds` serialization across operating systems
are not required. A future clean-machine validation with a committed
`renv.lock` will provide the strongest test of package-level reproducibility.

The three headline translation notebooks have already been tested locally
without their relevant caches. Their newly fitted non-XGBoost rows reproduced
the stored references exactly. This local result does not replace the final
fresh-clone test on a second machine.

## 10. GitHub and external data storage

Commit the source code, documentation, environment files, and small verified
pipeline exports to GitHub. In particular, keep the `.Rmd` and `.R` files,
`functions_all.R`, this guide, `README.md`, `IMPLEMENTATION.md`, the future
`renv.lock` and activation files, the six original-outcome XGBoost tuning
exports, the three XGBoost translation exports, the three combined translation
references, and the two Child four-unit exports. Include `ae98_sub.rds` only if
the underlying microdata terms permit redistribution of that record-level
derived file.

Do not commit knitr cache folders, `*_files/`, rendered HTML files, R session
state, or local `renv/library/` directories. Do not upload copies of published
papers, private guidelines, or other copyrighted reference documents as part
of the replication code.

The three raw `.dta` files may be committed only if their redistribution terms
permit it. If redistribution is uncertain or prohibited, place them in the
authorized external archive and give the reader exact acquisition instructions
and checksums. A general personal Dropbox link is not a substitute for stable
data provenance, but it can be used as a restricted delivery mechanism when
the repository documents the expected filenames and checksums.

## 11. Common errors

`Missing verified ... export` means that a required `.rds` file was not copied
with the repository. Restore that file or deliberately rebuild it using the
documented `TRUE` switch.

`uses tuning settings that differ` means that an export was produced with a
different seed, fold count, tuning budget, or tuning policy. Do not bypass this
check; restore the matching export or rebuild the analysis consistently.

An error after changing shared functions can be caused by a stale local knitr
cache. Preserve the current project first, restart R, and move aside only the
affected notebook cache. Never delete the verified `.rds` exports merely to
clear a cache.
