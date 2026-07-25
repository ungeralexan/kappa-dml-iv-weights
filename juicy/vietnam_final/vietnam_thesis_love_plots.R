#!/usr/bin/env Rscript

# Thesis-only balance figures for the Vietnam application.
# The presentation notebook remains unchanged. This script imports the saved
# tuned fits and creates two compact six-panel Love plots for the appendix.

suppressPackageStartupMessages({
  library(haven)
  library(cobalt)
  library(ggplot2)
})
options(warn = 1)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_arg)) {
  stop("Run this file with Rscript.")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
project_dir <- dirname(script_dir)
figure_dir <- file.path(project_dir, "pics")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(script_dir, "functions_all.R"))

# Reconstruct the common analysis sample and the two design matrices.
sipp <- read_dta(file.path(script_dir, "sipp.dta"))
sipp_clean <- sipp[
  sipp$rsncode != 999 &
    !is.na(sipp$educ) &
    !is.na(sipp$kwage),
]
sipp_clean$age <- sipp_clean$age_5
sipp_clean$age2 <- sipp_clean$age_5^2
sipp_clean$age3 <- sipp_clean$age_5^3
sipp_clean$age_fac <- factor(sipp_clean$age_5)

D <- as.integer(sipp_clean$nvstat)
Z <- as.integer(sipp_clean$rsncode)

X_cub <- model.matrix(~ age + age2 + age3, data = sipp_clean)
X_sat <- model.matrix(~ age_fac, data = sipp_clean)
X_dml_cub <- make_X(
  ~ 0 + age + I(age^2) + I(age^3),
  data = sipp_clean
)
X_dml_sat_full <- make_X(~ 0 + age_fac, data = sipp_clean)
X_dml_sat <- X_dml_sat_full[, -1, drop = FALSE]

stopifnot(
  nrow(sipp_clean) == 3027L,
  nrow(X_dml_cub) == length(D),
  nrow(X_dml_sat) == length(D)
)

# Import the cached tuned honest-forest outcome weights.
dml_cache <- Sys.glob(file.path(
  script_dir,
  "vietnam_presentation_4_cache",
  "html",
  "dml-fit-all_*.RData"
))
if (length(dml_cache) != 1L) {
  stop("Expected exactly one cached dml-fit-all object.")
}
dml_env <- new.env(parent = emptyenv())
invisible(lazyLoad(sub("\\.RData$", "", dml_cache), envir = dml_env))

# Import the tuned XGBoost outcome-weight reconstructions.
xgb_pliv <- readRDS(file.path(
  script_dir,
  "vietnam_xgb_pliv_tuning_export.rds"
))
xgb_wald <- readRDS(file.path(
  script_dir,
  "vietnam_xgb_tuning_export.rds"
))

# Recompute only the inexpensive kappa weights from the common sample.
kw_cub <- kappa_weights_bundle(Z, D, X_cub)
kw_sat <- kappa_weights_bundle(Z, D, X_sat)

weights_cub <- list(
  "PLR-IV\nTuned forest" = dml_env$w_cub_tuned$plriv,
  "Wald-AIPW\nTuned forest" = dml_env$w_cub_tuned$waldaipw,
  "PLR-IV\nTuned XGBoost" = xgb_pliv$omega_vectors[["Cubic age"]]$tuned,
  "Wald-AIPW\nTuned XGBoost" = xgb_wald$omega_vectors[["Cubic age"]]$tuned,
  "Kappa\nCBPS" = kw_cub$tau_cb_u,
  "Kappa\nMLE" = kw_cub$tau_ml_u
)

weights_sat <- list(
  "PLR-IV\nTuned forest" = dml_env$w_sat_tuned$plriv,
  "Wald-AIPW\nTuned forest" = dml_env$w_sat_tuned$waldaipw,
  "PLR-IV\nTuned XGBoost" = xgb_pliv$omega_vectors[["Saturated age"]]$tuned,
  "Wald-AIPW\nTuned XGBoost" = xgb_wald$omega_vectors[["Saturated age"]]$tuned,
  "Kappa\nCBPS" = kw_sat$tau_cb_u,
  "Kappa\nMLE" = kw_sat$tau_ml_u
)

# Extract the same balance statistic used by cobalt::love.plot().
balance_frame <- function(label, w, X, variable_labels) {
  bt <- cobalt::bal.tab(
    x = X,
    treat = D,
    weights = w * (2 * D - 1),
    continuous = "std",
    binary = "raw",
    s.d.denom = "pooled",
    un = TRUE,
    quick = FALSE
  )
  balance <- bt$Balance
  stopifnot(
    all(c("Diff.Un", "Diff.Adj") %in% names(balance)),
    nrow(balance) == length(variable_labels)
  )
  data.frame(
    Estimator = label,
    Variable = rep(variable_labels, 2L),
    Sample = rep(c("Unadjusted", "Adjusted"), each = nrow(balance)),
    ASMD = abs(c(balance$Diff.Un, balance$Diff.Adj)),
    stringsAsFactors = FALSE
  )
}

collect_balance <- function(weight_list, X, variable_labels) {
  labels <- names(weight_list)
  out <- do.call(
    rbind,
    lapply(seq_along(weight_list), function(i) {
      balance_frame(labels[[i]], weight_list[[i]], X, variable_labels)
    })
  )
  out$Estimator <- factor(out$Estimator, levels = labels)
  out$Variable <- factor(
    out$Variable,
    levels = rev(variable_labels)
  )
  out$Sample <- factor(
    out$Sample,
    levels = c("Unadjusted", "Adjusted")
  )
  rownames(out) <- NULL
  out
}

cubic_labels <- c("Age", "Age^2", "Age^3")
saturated_labels <- sub("^age_fac", "Age ", colnames(X_dml_sat))

balance_cub <- collect_balance(weights_cub, X_dml_cub, cubic_labels)
balance_sat <- collect_balance(weights_sat, X_dml_sat, saturated_labels)

make_thesis_love <- function(df, title, x_label, x_limit) {
  ggplot(
    df,
    aes(
      x = ASMD,
      y = Variable,
      colour = Sample,
      shape = Sample,
      group = Sample
    )
  ) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.45) +
    geom_vline(
      xintercept = 0.1,
      colour = "grey30",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.8) +
    facet_wrap(~ Estimator, ncol = 3) +
    scale_colour_manual(
      values = c("Unadjusted" = "#440154", "Adjusted" = "#FDE725")
    ) +
    scale_shape_manual(
      values = c("Unadjusted" = 16, "Adjusted" = 17)
    ) +
    scale_x_continuous(
      limits = c(0, x_limit),
      expand = expansion(mult = c(0.015, 0.03))
    ) +
    labs(
      title = title,
      x = x_label,
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11,
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(1.1, "lines"),
      axis.text = element_text(colour = "black"),
      legend.position = "bottom",
      legend.box.margin = margin(t = 4),
      plot.margin = margin(8, 12, 8, 8)
    )
}

plot_cub <- make_thesis_love(
  balance_cub,
  "Covariate balance under the cubic age specification",
  "Absolute standardized mean difference",
  x_limit = 0.58
)
plot_sat <- make_thesis_love(
  balance_sat,
  "Covariate balance under the saturated age specification",
  "Absolute mean difference",
  x_limit = 0.105
)

save_plot <- function(plot, stem, width, height) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = pdf,
    useDingbats = FALSE
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

save_plot(
  plot_cub,
  "vietnam_love_cubic_thesis",
  width = 7,
  height = 5.1
)
save_plot(
  plot_sat,
  "vietnam_love_saturated_thesis",
  width = 7,
  height = 6.6
)

balance_summary <- rbind(
  transform(
    subset(balance_cub, Sample == "Adjusted"),
    Specification = "Cubic age"
  ),
  transform(
    subset(balance_sat, Sample == "Adjusted"),
    Specification = "Saturated age"
  )
)
balance_summary <- aggregate(
  ASMD ~ Specification + Estimator,
  data = balance_summary,
  FUN = max
)
names(balance_summary)[names(balance_summary) == "ASMD"] <- "Max_adjusted_ASMD"
write.csv(
  balance_summary,
  file.path(figure_dir, "vietnam_love_balance_summary.csv"),
  row.names = FALSE
)

message("Created thesis Love plots and balance summary in: ", figure_dir)
