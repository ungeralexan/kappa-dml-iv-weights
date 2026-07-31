# Focused covariate-balance figures for the Card application.
#
# This script is sourced at the end of card_test.Rmd, after candidate_cells and
# the four design objects have been constructed. It creates two compact PNG
# figures in ../pics and deliberately excludes outcome-weight vectors that do
# not pass the estimator-reproduction check.

card_balance_one_covariate <- function(x, D, oriented_weights = NULL) {
  x <- as.numeric(x)
  D <- as.integer(D)

  is_binary <- all(stats::na.omit(unique(x)) %in% c(0, 1))

  if (is.null(oriented_weights)) {
    mean_t <- mean(x[D == 1], na.rm = TRUE)
    mean_c <- mean(x[D == 0], na.rm = TRUE)
  } else {
    a <- as.numeric(oriented_weights)
    mean_t <- sum(a[D == 1] * x[D == 1], na.rm = TRUE) /
      sum(a[D == 1], na.rm = TRUE)
    mean_c <- sum(a[D == 0] * x[D == 0], na.rm = TRUE) /
      sum(a[D == 0], na.rm = TRUE)
  }

  difference <- abs(mean_t - mean_c)
  if (!is_binary) {
    pooled_sd <- sqrt(
      (stats::var(x[D == 1], na.rm = TRUE) +
         stats::var(x[D == 0], na.rm = TRUE)) / 2
    )
    difference <- difference / pooled_sd
  }

  c(balance = difference, binary = as.numeric(is_binary))
}

card_balance_variable_labels <- function(variable_names, binary_flags) {
  clean <- c(
    exper = "Experience",
    expersq = "Experience squared",
    black = "Black",
    south = "South",
    smsa = "SMSA",
    smsa66 = "SMSA, 1966",
    south66 = "South, 1966",
    reg661 = "Region 1, 1966",
    reg662 = "Region 2, 1966",
    reg663 = "Region 3, 1966",
    reg664 = "Region 4, 1966",
    reg665 = "Region 5, 1966",
    reg666 = "Region 6, 1966",
    reg667 = "Region 7, 1966",
    reg668 = "Region 8, 1966"
  )

  labels <- unname(clean[variable_names])
  labels[is.na(labels)] <- variable_names[is.na(labels)]
  ifelse(binary_flags, paste0(labels, "*"), labels)
}

card_selected_balance_data <- function(cell_name, treatment_label, D, X,
                                       candidate_cells) {
  selected <- c(
    "CBPS kappa" = "τ_cb_u (kappa, CBPS)",
    "MLE kappa" = "τ_ml_u (kappa, MLE)",
    "PLR-IV forest" = "PLR-IV, tuned forest",
    "Wald-AIPW forest" = "Wald-AIPW, tuned forest"
  )

  X <- as.matrix(X)
  variables <- colnames(X)
  if (is.null(variables)) variables <- paste0("X", seq_len(ncol(X)))

  raw_stats <- t(vapply(
    seq_len(ncol(X)),
    function(j) card_balance_one_covariate(X[, j], D),
    numeric(2)
  ))

  out <- data.frame(
    Treatment = treatment_label,
    Variable = variables,
    Series = "Unadjusted",
    Balance = raw_stats[, "balance"],
    Binary = as.logical(raw_stats[, "binary"]),
    stringsAsFactors = FALSE
  )

  candidates <- candidate_cells[[cell_name]]$candidates
  for (display_name in names(selected)) {
    candidate_name <- unname(selected[[display_name]])
    candidate <- candidates[[candidate_name]]
    if (is.null(candidate) || !isTRUE(candidate$ok)) {
      stop("Eligible balance weight not found: ", cell_name, " / ", candidate_name)
    }

    oriented <- as.numeric(candidate$w) * (2 * D - 1)
    adjusted_stats <- t(vapply(
      seq_len(ncol(X)),
      function(j) card_balance_one_covariate(X[, j], D, oriented),
      numeric(2)
    ))

    out <- rbind(
      out,
      data.frame(
        Treatment = treatment_label,
        Variable = variables,
        Series = display_name,
        Balance = adjusted_stats[, "balance"],
        Binary = as.logical(adjusted_stats[, "binary"]),
        stringsAsFactors = FALSE
      )
    )
  }

  binary_map <- stats::setNames(out$Binary[out$Series == "Unadjusted"],
                                out$Variable[out$Series == "Unadjusted"])
  out$Variable_label <- card_balance_variable_labels(
    out$Variable,
    unname(binary_map[out$Variable])
  )
  out
}

card_plot_balance <- function(balance_data, specification, output_file,
                              width, height) {
  series_order <- c(
    "Unadjusted", "CBPS kappa", "MLE kappa",
    "PLR-IV forest", "Wald-AIPW forest"
  )
  offsets <- c(
    "Unadjusted" = -0.24,
    "CBPS kappa" = -0.12,
    "MLE kappa" = 0,
    "PLR-IV forest" = 0.12,
    "Wald-AIPW forest" = 0.24
  )

  variable_order <- unique(balance_data$Variable_label)
  variable_order <- rev(variable_order)
  variable_index <- stats::setNames(seq_along(variable_order), variable_order)

  balance_data$Series <- factor(balance_data$Series, levels = series_order)
  balance_data$Treatment <- factor(
    balance_data$Treatment,
    levels = c("Schooling beyond high school", "At least 16 years")
  )
  balance_data$Variable_label <- factor(
    balance_data$Variable_label,
    levels = variable_order
  )
  balance_data$y <- unname(variable_index[as.character(balance_data$Variable_label)]) +
    unname(offsets[as.character(balance_data$Series)])

  colors <- c(
    "Unadjusted" = "#666666",
    "CBPS kappa" = "#000000",
    "MLE kappa" = "#E69F00",
    "PLR-IV forest" = "#0072B2",
    "Wald-AIPW forest" = "#009E73"
  )
  shapes <- c(
    "Unadjusted" = 1,
    "CBPS kappa" = 18,
    "MLE kappa" = 17,
    "PLR-IV forest" = 15,
    "Wald-AIPW forest" = 16
  )

  plot <- ggplot2::ggplot(
    balance_data,
    ggplot2::aes(x = Balance, y = y, color = Series, shape = Series)
  ) +
    ggplot2::geom_hline(
      yintercept = seq_along(variable_order),
      color = "grey92",
      linewidth = 0.35
    ) +
    ggplot2::geom_vline(
      xintercept = 0.10,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +
    ggplot2::geom_point(size = 2.2, stroke = 0.8) +
    ggplot2::facet_wrap(~Treatment, nrow = 1) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(variable_order),
      labels = variable_order,
      expand = ggplot2::expansion(add = 0.55)
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.05))
    ) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shapes, drop = FALSE) +
    ggplot2::labs(
      title = paste0("Outcome-weight covariate balance: ", specification,
                     " conditioning set"),
      subtitle = "Selected exact outcome-weight representations",
      x = "Absolute covariate mean difference",
      y = NULL,
      color = NULL,
      shape = NULL,
      caption = paste(
        "Dashed line: 0.10. * Binary covariates use raw differences;",
        "continuous covariates use pooled-SD standardized differences."
      )
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11.5),
      plot.subtitle = ggplot2::element_text(size = 9.5, color = "grey35"),
      plot.caption = ggplot2::element_text(size = 8, color = "grey35", hjust = 0),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8.5, color = "black"),
      axis.text.x = ggplot2::element_text(size = 8.5, color = "black"),
      axis.title.x = ggplot2::element_text(size = 9.5),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8.5),
      legend.key.width = grid::unit(1.2, "lines"),
      panel.spacing = grid::unit(1.0, "lines"),
      plot.margin = ggplot2::margin(8, 12, 6, 8)
    )

  ggplot2::ggsave(
    filename = output_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 320,
    bg = "white"
  )

  invisible(plot)
}

create_card_balance_figures <- function(candidate_cells, D1, D2,
                                        X_dml_card, X_dml_kit,
                                        output_dir = "../pics") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  card_data <- rbind(
    card_selected_balance_data(
      "Card / some college", "Schooling beyond high school",
      D1, X_dml_card, candidate_cells
    ),
    card_selected_balance_data(
      "Card / college completion", "At least 16 years",
      D2, X_dml_card, candidate_cells
    )
  )

  kitagawa_data <- rbind(
    card_selected_balance_data(
      "Kitagawa / some college", "Schooling beyond high school",
      D1, X_dml_kit, candidate_cells
    ),
    card_selected_balance_data(
      "Kitagawa / college completion", "At least 16 years",
      D2, X_dml_kit, candidate_cells
    )
  )

  card_plot_balance(
    card_data,
    "Card",
    file.path(output_dir, "card_balance_card.png"),
    width = 8.4,
    height = 5.7
  )
  card_plot_balance(
    kitagawa_data,
    "Kitagawa",
    file.path(output_dir, "card_balance_kitagawa.png"),
    width = 8.4,
    height = 3.8
  )

  invisible(list(card = card_data, kitagawa = kitagawa_data))
}
