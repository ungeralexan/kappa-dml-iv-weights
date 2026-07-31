## Export the already-computed Card figure used by the current thesis.
## This script does not fit or tune any learner.
## Run from the project root.

source_dir <- normalizePath("card_final", mustWork = TRUE)
target_dir <- normalizePath("pics", mustWork = TRUE)

figure_map <- c(
  "card_descriptives_1_files/figure-html/ps-dist-plot-1.png" =
    "card_propensity_support.png"
)

for (relative_source in names(figure_map)) {
  copied <- file.copy(
    file.path(source_dir, relative_source),
    file.path(target_dir, unname(figure_map[[relative_source]])),
    overwrite = TRUE
  )
  if (!copied) stop("Could not copy ", relative_source)
}

message("Exported ", length(figure_map), " current Card thesis figure to ", target_dir)
