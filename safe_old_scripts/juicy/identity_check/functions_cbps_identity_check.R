# Focused finite-sample identity diagnostics for the kappa estimators.
# The estimator and propensity-score implementations remain in functions_all.R;
# this file only assembles them into the comparison used in the thesis.

require_kappa_identity_dependencies <- function() {
  required <- c(
    "logit_mle", "cbps", "tau_u", "tau_a10", "tau_unnorm",
    "kappa_weights"
  )
  missing <- required[!vapply(required, exists, logical(1), mode = "function")]
  if (length(missing)) {
    stop(
      "Source functions_all.R before functions_cbps_identity_check.R. Missing: ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}

has_constant_column <- function(X, tolerance = 1e-12) {
  X <- as.matrix(X)
  any(apply(X, 2, function(x) all(is.finite(x)) && max(abs(x - 1)) <= tolerance))
}

kappa_estimators_at_p <- function(Y, Z, D, p) {
  c(
    tau_u   = tau_u(Y, Z, D, p),
    tau_a10 = tau_a10(Y, Z, D, p),
    tau_a   = tau_unnorm(Y, Z, D, p, "a"),
    tau_a1  = tau_unnorm(Y, Z, D, p, "a1"),
    tau_a0  = tau_unnorm(Y, Z, D, p, "a0")
  )
}

kappa_denominator_diagnostics <- function(Z, D, p) {
  kw <- kappa_weights(Z, D, p)
  r <- (Z - p) / (p * (1 - p))
  c(
    instrument_mass_Z1 = sum(Z / p),
    instrument_mass_Z0 = sum((1 - Z) / (1 - p)),
    instrument_mass_difference = sum(Z / p) - sum((1 - Z) / (1 - p)),
    intercept_balancing_sum = sum(r),
    sum_kappa1 = sum(kw$kappa1),
    sum_kappa0 = sum(kw$kappa0),
    kappa1_minus_kappa0 = sum(kw$kappa1) - sum(kw$kappa0),
    sum_kappa = sum(kw$kappa)
  )
}

pairwise_minimum_gap <- function(x) {
  gaps <- abs(outer(as.numeric(x), as.numeric(x), "-"))
  gaps[lower.tri(gaps)] |> min()
}

check_cbps_finite_sample_identity <- function(Y, Z, D, X,
                                              identity_tolerance = 1e-7,
                                              cbps_tolerance = 1e-9,
                                              cbps_max_iter = 5000) {
  require_kappa_identity_dependencies()

  Y <- as.numeric(Y)
  Z <- as.numeric(Z)
  D <- as.numeric(D)
  X <- as.matrix(X)
  n <- length(Y)

  if (!all(c(length(Z), length(D), nrow(X)) == n)) {
    stop("Y, Z, D, and X must contain the same number of observations.")
  }
  if (!has_constant_column(X)) {
    stop("The CBPS identity requires an explicit constant column in X.")
  }

  p_ml <- logit_mle(Z, X)
  cb_fit <- cbps(Z, X, tol = cbps_tolerance, max_iter = cbps_max_iter)
  p_cb <- as.numeric(cb_fit$p)

  estimates_ml <- kappa_estimators_at_p(Y, Z, D, p_ml)
  estimates_cb <- kappa_estimators_at_p(Y, Z, D, p_cb)

  # tau_t is exactly the same estimator as tau_a1, so it is not duplicated here.
  identity_names <- c("tau_u", "tau_a10", "tau_a1", "tau_a0")
  cb_reference <- unname(estimates_cb["tau_u"])
  cb_gaps <- estimates_cb[identity_names] - cb_reference
  cb_max_gap <- max(abs(cb_gaps))

  estimate_table <- data.frame(
    Estimator = c("tau_u", "tau_a10", "tau_a", "tau_a1 (= tau_t)", "tau_a0"),
    MLE = unname(estimates_ml[c("tau_u", "tau_a10", "tau_a", "tau_a1", "tau_a0")]),
    CBPS = unname(estimates_cb[c("tau_u", "tau_a10", "tau_a", "tau_a1", "tau_a0")]),
    In_CBPS_equivalence_class = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    CBPS_gap_from_tau_u = c(
      0,
      unname(estimates_cb["tau_a10"] - cb_reference),
      unname(estimates_cb["tau_a"] - cb_reference),
      unname(estimates_cb["tau_a1"] - cb_reference),
      unname(estimates_cb["tau_a0"] - cb_reference)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  list(
    estimates = estimate_table,
    propensity = list(ml = p_ml, cb = p_cb),
    denominator_diagnostics = rbind(
      MLE = kappa_denominator_diagnostics(Z, D, p_ml),
      CBPS = kappa_denominator_diagnostics(Z, D, p_cb)
    ),
    checks = list(
      constant_in_X = TRUE,
      cbps_converged = isTRUE(cb_fit$converged),
      cbps_max_moment = as.numeric(cb_fit$max_moment),
      cbps_identity_tolerance = identity_tolerance,
      cbps_max_equivalence_gap = cb_max_gap,
      cbps_identity_pass = is.finite(cb_max_gap) && cb_max_gap <= identity_tolerance,
      cbps_tau_a_gap = unname(estimates_cb["tau_a"] - cb_reference),
      ml_minimum_pairwise_gap = pairwise_minimum_gap(estimates_ml),
      ml_all_five_numerically_distinct =
        pairwise_minimum_gap(estimates_ml) > identity_tolerance
    )
  )
}
