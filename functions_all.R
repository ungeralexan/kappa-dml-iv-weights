# Convert predictors to probabilities.
safe_logit <- function(eta) {
  eta <- pmin(pmax(eta, -35), 35)
  1 / (1 + exp(-eta))
}

# Standardize the design matrix.
prep_design_for_mest <- function(X) {
  X <- as.matrix(X)
  is_intercept <- apply(X, 2, function(v) all(abs(v - 1) < 1e-12))
  X_new <- X
  for (j in seq_len(ncol(X))) {
    if (!is_intercept[j]) {
      mu  <- mean(X[, j])
      sdj <- sd(X[, j])
      if (is.finite(sdj) && sdj > 1e-12)
        X_new[, j] <- (X[, j] - mu) / sdj
    }
  }
  X_new
}

# Fit MLE instrument propensities.
logit_mle <- function(Z, X) {
  df  <- data.frame(Z = Z, X)
  fit <- glm(Z ~ . - 1, data = df, family = binomial(link = "logit"))
  fitted.values(fit)
}

# Prepare MLE inputs for standard errors.
fit_logit_alpha <- function(Z, X) {
  Z     <- as.numeric(Z)
  X     <- prep_design_for_mest(X)
  fit   <- glm.fit(x = X, y = Z, family = binomial(link = "logit"))
  alpha <- as.numeric(coef(fit));  alpha[is.na(alpha)] <- 0
  p     <- as.vector(pmin(pmax(safe_logit(X %*% alpha), 1e-8), 1 - 1e-8))
  list(alpha = alpha, p = p, X_used = X)
}

# Fit balancing instrument propensities.
fit_cbps_alpha <- function(Z, X, tol = 1e-9, max_iter = 5000) {
  Z <- as.numeric(Z)
  X <- prep_design_for_mest(X)
  n <- length(Z)
  k <- ncol(X)

  b <- tryCatch(
    as.numeric(glm.fit(x = X, y = Z, family = binomial())$coefficients),
    error = function(e) rep(0, k)
  )
  b[is.na(b)] <- 0

  # Calculate CBPS moments.
  moment_fn <- function(b) {
    p <- as.vector(pmin(pmax(safe_logit(X %*% b), 1e-8), 1 - 1e-8))
    colMeans(as.vector((Z - p) / (p * (1 - p))) * X)
  }

  # Calculate the CBPS Jacobian.
  jac_fn <- function(b) {
    p <- as.vector(pmin(pmax(safe_logit(X %*% b), 1e-8), 1 - 1e-8))
    w <- as.vector(-Z * (1 - p) / p - (1 - Z) * p / (1 - p))
    crossprod(X, w * X) / n
  }

  best_b    <- b
  best_norm <- max(abs(moment_fn(b)))
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    m <- moment_fn(b);  m_norm <- max(abs(m))
    if (m_norm < best_norm) { best_norm <- m_norm;  best_b <- b }
    if (m_norm < tol)       { converged <- TRUE;    best_b <- b;  break }

    J    <- jac_fn(b)
    step <- tryCatch(qr.solve(J, -m), error = function(e) NULL)

    if (is.null(step) || any(!is.finite(step))) {
      for (ridge in c(1e-10, 1e-8, 1e-6, 1e-4)) {
        step <- tryCatch(solve(J + ridge * diag(k), -m), error = function(e) NULL)
        if (!is.null(step) && all(is.finite(step))) break
      }
    }
    if (is.null(step) || any(!is.finite(step))) break

    alpha_step <- 1
    for (j in seq_len(50)) {
      b_new <- b + alpha_step * step
      if (is.finite(max(abs(moment_fn(b_new)))) &&
          max(abs(moment_fn(b_new))) < m_norm) { b <- b_new;  break }
      alpha_step <- alpha_step * 0.5
    }
  }

  b <- best_b
  p <- as.vector(pmin(pmax(safe_logit(X %*% b), 1e-8), 1 - 1e-8))
  list(alpha = as.numeric(b), p = p, X_used = X,
       converged = converged, max_moment = best_norm)
}

# Return the complete CBPS fit.
cbps <- function(Z, X, tol = 1e-9, max_iter = 5000, verbose = FALSE) {
  fit_cbps_alpha(Z, X, tol = tol, max_iter = max_iter)
}

# Extract CBPS propensities.
get_cbps_p <- function(Z, X) {
  out <- cbps(Z, X)
  if (is.list(out) && !is.null(out$p)) return(as.vector(out$p))
  as.vector(out)
}

# Calculate the three kappa weights.
kappa_weights <- function(Z, D, p) {
  list(
    kappa  = 1 - D * (1 - Z) / (1 - p) - (1 - D) * Z / p,
    kappa1 = D * (Z - p) / (p * (1 - p)),
    kappa0 = (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
  )
}

# Calculate the Uysal estimator.
tau_u <- function(Y, Z, D, p) {
  s1 <- sum(Z / p)
  s0 <- sum((1 - Z) / (1 - p))
  numerator   <- sum(Y * Z / p) / s1 - sum(Y * (1 - Z) / (1 - p)) / s0
  denominator <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  numerator / denominator
}

# Calculate the Abadie-Cattaneo estimator.
tau_a10 <- function(Y, Z, D, p) {
  kw <- kappa_weights(Z, D, p)
  sum(kw$kappa1 * Y) / sum(kw$kappa1) - sum(kw$kappa0 * Y) / sum(kw$kappa0)
}

# Calculate an unnormalized kappa estimator.
tau_unnorm <- function(Y, Z, D, p, which = "a") {
  kw        <- kappa_weights(Z, D, p)
  numerator <- mean(Y * (Z - p) / (p * (1 - p)))
  denom_val <- switch(which,
                      "a"  = mean(kw$kappa),
                      "a1" = mean(kw$kappa1),
                      "a0" = mean(kw$kappa0)
  )
  numerator / denom_val
}

# Construct kappa outcome weights.
kappa_outcome_weights <- function(Z, D, p) {
  n  <- length(Z)
  kw <- kappa_weights(Z, D, p)

  s1  <- sum(Z / p)
  s0  <- sum((1 - Z) / (1 - p))
  dD  <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  w_u <- (Z / p / s1 - (1 - Z) / (1 - p) / s0) / dD

  w_a10 <- kw$kappa1 / sum(kw$kappa1) - kw$kappa0 / sum(kw$kappa0)

  num_w <- (Z - p) / (p * (1 - p)) / n

  list(
    w_u   = as.vector(w_u),
    w_a10 = as.vector(w_a10),
    w_a   = as.vector(num_w / mean(kw$kappa)),
    w_a1  = as.vector(num_w / mean(kw$kappa1)),
    w_a0  = as.vector(num_w / mean(kw$kappa0))
  )
}

# Approximate a numerical Jacobian.
num_jacobian <- function(f, theta, eps = 1e-6) {
  theta <- as.numeric(theta)
  f0 <- f(theta)
  m  <- length(f0)
  k  <- length(theta)
  J  <- matrix(NA_real_, nrow = m, ncol = k)

  for (j in seq_len(k)) {
    h  <- eps * (abs(theta[j]) + 1)
    tp <- theta;  tp[j] <- tp[j] + h
    tm <- theta;  tm[j] <- tm[j] - h
    J[, j] <- (f(tp) - f(tm)) / (2 * h)
  }
  J
}

# Invert a matrix safely.
matrix_inverse_safe <- function(A, tol = 1e-10) {
  inv <- tryCatch(solve(A), error = function(e) NULL)
  if (!is.null(inv) && all(is.finite(inv))) return(inv)

  inv <- tryCatch(qr.solve(A), error = function(e) NULL)
  if (!is.null(inv) && all(is.finite(inv))) return(inv)

  for (ridge in c(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)) {
    inv <- tryCatch(solve(A + ridge * diag(ncol(A))), error = function(e) NULL)
    if (!is.null(inv) && all(is.finite(inv))) return(inv)
  }

  sv    <- svd(A)
  d     <- sv$d
  d_inv <- ifelse(d > tol * max(d), 1 / d, 0)
  sv$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(sv$u)
}

# Calculate a sandwich standard error.
sandwich_se_mest <- function(moment_matrix_fn, theta_hat, tau_index) {
  theta_hat    <- as.numeric(theta_hat)
  psi_hat      <- moment_matrix_fn(theta_hat)
  n            <- nrow(psi_hat)

  A            <- num_jacobian(function(th) colMeans(moment_matrix_fn(th)), theta_hat)
  psi_centered <- scale(psi_hat, center = TRUE, scale = FALSE)
  V            <- crossprod(psi_centered) / n

  A_inv        <- matrix_inverse_safe(A)
  vcov_theta   <- A_inv %*% V %*% t(A_inv) / n

  se2 <- vcov_theta[tau_index, tau_index]
  if (!is.finite(se2)) return(NA_real_)
  sqrt(abs(se2))
}

# Construct propensity-score moments.
alpha_moment_matrix <- function(Z, p, X_used, method) {
  p <- as.vector(p);  Z <- as.vector(Z)
  if (method == "ml") return(as.vector(Z - p) * X_used)
  if (method == "cb") return(as.vector((Z - p) / (p * (1 - p))) * X_used)
  stop("method must be 'ml' or 'cb'")
}

# Calculate one analytical kappa standard error.
kappa_analytic_se_one <- function(Y, Z, D, X, estimator, method = "ml") {
  Y <- as.numeric(Y);  Z <- as.numeric(Z)
  D <- as.numeric(D);  X <- as.matrix(X)
  n <- length(Y)

  fit       <- if (method == "ml") fit_logit_alpha(Z, X) else fit_cbps_alpha(Z, X)
  alpha_hat <- fit$alpha
  X_used    <- fit$X_used
  k         <- length(alpha_hat)

  # Select the propensity moments.
  am <- function(p) alpha_moment_matrix(Z, p, X_used, method)

  if (estimator == "u") {
    p   <- fit$p
    s1  <- sum(Z / p);  s0 <- sum((1 - Z) / (1 - p))
    mu1 <- sum(Z * Y / p) / s1;            mu0 <- sum((1 - Z) * Y / (1 - p)) / s0
    m1  <- sum(Z * D / p) / s1;            m0  <- sum((1 - Z) * D / (1 - p)) / s0
    tau <- (mu1 - mu0) / (m1 - m0)

    theta_hat <- c(alpha_hat, mu1, mu0, m1, m0, tau)

    # Stack the Uysal moments.
    moment_fn <- function(theta) {
      a   <- theta[seq_len(k)]
      mu1 <- theta[k+1]; mu0 <- theta[k+2]
      m1  <- theta[k+3]; m0  <- theta[k+4]
      tau <- theta[k+5]
      p   <- as.vector(pmin(pmax(safe_logit(X_used %*% a), 1e-8), 1-1e-8))
      cbind(am(p),
            psi_mu1 = Z * (Y - mu1) / p,
            psi_mu0 = (1 - Z) * (Y - mu0) / (1 - p),
            psi_m1  = Z * (D - m1) / p,
            psi_m0  = (1 - Z) * (D - m0) / (1 - p),
            psi_tau = (mu1 - mu0) / (m1 - m0) - tau)
    }
    return(sandwich_se_mest(moment_fn, theta_hat, length(theta_hat)))
  }

  if (estimator %in% c("a", "a1", "a0")) {
    p     <- fit$p
    kw    <- kappa_weights(Z, D, p)
    Delta <- mean(Y * (Z - p) / (p * (1 - p)))
    Gamma <- switch(estimator,
                    "a"  = mean(kw$kappa),
                    "a1" = mean(kw$kappa1),
                    "a0" = mean(kw$kappa0))
    tau       <- Delta / Gamma
    theta_hat <- c(alpha_hat, Delta, Gamma, tau)

    # Stack the unnormalized moments.
    moment_fn <- function(theta) {
      a     <- theta[seq_len(k)]
      Delta <- theta[k+1]; Gamma <- theta[k+2]; tau <- theta[k+3]
      p     <- as.vector(pmin(pmax(safe_logit(X_used %*% a), 1e-8), 1-1e-8))
      psi_Delta <- Z * Y / p - (1 - Z) * Y / (1 - p) - Delta
      psi_Gamma <- switch(estimator,
                          "a"  = 1 - (1-Z)*D/(1-p) - Z*(1-D)/p - Gamma,
                          "a1" = Z*D/p - (1-Z)*D/(1-p) - Gamma,
                          "a0" = Z*(D-1)/p - (1-Z)*(D-1)/(1-p) - Gamma)
      cbind(am(p),
            psi_Delta = psi_Delta,
            psi_Gamma = psi_Gamma,
            psi_tau   = Delta / Gamma - tau)
    }
    return(sandwich_se_mest(moment_fn, theta_hat, length(theta_hat)))
  }

  if (estimator == "a10") {
    p   <- fit$p
    kw  <- kappa_weights(Z, D, p)
    Delta1 <- mean(kw$kappa1 * Y);  Gamma1 <- mean(kw$kappa1)
    Delta0 <- mean(kw$kappa0 * Y);  Gamma0 <- mean(kw$kappa0)
    tau    <- Delta1/Gamma1 - Delta0/Gamma0
    theta_hat <- c(alpha_hat, Delta1, Gamma1, Delta0, Gamma0, tau)

    # Stack the Abadie-Cattaneo moments.
    moment_fn <- function(theta) {
      a      <- theta[seq_len(k)]
      Delta1 <- theta[k+1]; Gamma1 <- theta[k+2]
      Delta0 <- theta[k+3]; Gamma0 <- theta[k+4]; tau <- theta[k+5]
      p      <- as.vector(pmin(pmax(safe_logit(X_used %*% a), 1e-8), 1-1e-8))
      k1_i   <- D * (Z - p) / (p * (1 - p))
      k0_i   <- (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
      cbind(am(p),
            psi_Delta1 = k1_i * Y - Delta1,
            psi_Gamma1 = Z*D/p - (1-Z)*D/(1-p) - Gamma1,
            psi_Delta0 = k0_i * Y - Delta0,
            psi_Gamma0 = Z*(D-1)/p - (1-Z)*(D-1)/(1-p) - Gamma0,
            psi_tau    = Delta1/Gamma1 - Delta0/Gamma0 - tau)
    }
    return(sandwich_se_mest(moment_fn, theta_hat, length(theta_hat)))
  }

  stop("Unknown estimator: must be one of 'u', 'a10', 'a', 'a1', 'a0'")
}

# Catch failed kappa standard errors.
safe_kappa_se <- function(Y, Z, D, X_mat, estimator, method) {
  tryCatch(
    kappa_analytic_se_one(Y, Z, D, X_mat, estimator = estimator, method = method),
    error = function(e) {
      warning(sprintf("SE failed: estimator=%s, method=%s. %s",
                      estimator, method, e$message))
      NA_real_
    }
  )
}

# Calculate all kappa point estimates.
kappa_point_estimates <- function(Y, Z, D, X_mat) {
  Y <- as.numeric(Y); Z <- as.numeric(Z)
  D <- as.numeric(D); X_mat <- as.matrix(X_mat)

  p_ml <- logit_mle(Z, X_mat)
  p_cb <- get_cbps_p(Z, X_mat)

  c(
    tau_cb_u   = tau_u(Y, Z, D, p_cb),
    tau_ml_u   = tau_u(Y, Z, D, p_ml),
    tau_ml_a10 = tau_a10(Y, Z, D, p_ml),
    tau_ml_a   = tau_unnorm(Y, Z, D, p_ml, "a"),
    tau_ml_t   = tau_unnorm(Y, Z, D, p_ml, "a1"),
    tau_ml_a0  = tau_unnorm(Y, Z, D, p_ml, "a0")
  )
}

# Calculate all kappa estimates and standard errors.
kappa_analytic_se_all <- function(Y, Z, D, X_mat) {
  Y <- as.numeric(Y); Z <- as.numeric(Z)
  D <- as.numeric(D); X_mat <- as.matrix(X_mat)

  estimates <- kappa_point_estimates(Y, Z, D, X_mat)

  se <- c(
    tau_cb_u   = safe_kappa_se(Y, Z, D, X_mat, "u",   "cb"),
    tau_ml_u   = safe_kappa_se(Y, Z, D, X_mat, "u",   "ml"),
    tau_ml_a10 = safe_kappa_se(Y, Z, D, X_mat, "a10", "ml"),
    tau_ml_a   = safe_kappa_se(Y, Z, D, X_mat, "a",   "ml"),
    tau_ml_t   = safe_kappa_se(Y, Z, D, X_mat, "a1",  "ml"),
    tau_ml_a0  = safe_kappa_se(Y, Z, D, X_mat, "a0",  "ml")
  )

  list(estimates = estimates, se = se)
}

# Fit the 2SLS benchmark.
run_2sls <- function(Y, D, Z, X_df, endog_name = "D") {
  cov_names <- names(X_df)[names(X_df) != "(Intercept)"]
  df <- data.frame(Y = Y, D = D, Z = Z, X_df)
  names(df)[names(df) == "D"] <- endog_name
  names(df)[names(df) == "Z"] <- "instrument"

  fml <- if (length(cov_names) == 0) {
    as.formula(paste("Y ~", endog_name, "| instrument"))
  } else {
    cs <- paste(cov_names, collapse = " + ")
    as.formula(paste("Y ~", endog_name, "+", cs, "| instrument +", cs))
  }

  fit    <- ivreg(fml, data = df)
  vcov_r <- vcovHC(fit, type = "HC1")
  ct     <- coeftest(fit, vcov = vcov_r)
  list(
    coef = ct[endog_name, "Estimate"],
    se   = ct[endog_name, "Std. Error"],
    pval = ct[endog_name, "Pr(>|t|)"]
  )
}

# Calculate weight summary statistics.
weight_statistics <- function(w) {
  w <- as.numeric(w)
  if (!length(w)) stop("Weight vector must contain at least one observation.")

  n <- length(w)
  if (any(!is.finite(w))) {
    return(c(
      n = n, Min = NA_real_, Max = NA_real_, Neg_share = NA_real_,
      Pct_neg = NA_real_, Sum_top10 = NA_real_, Sum_w = NA_real_,
      Sum_abs_w = NA_real_, Max_abs_w = NA_real_, ESS_kish = NA_real_,
      ESS_mod = NA_real_
    ))
  }

  sw2 <- sum(w^2)
  ord <- sort(w)
  p90 <- floor(0.9 * n)
  sum_top10 <- if (p90 < n) sum(ord[(p90 + 1):n]) else sum(ord)

  c(
    n          = n,
    Min        = min(w),
    Max        = max(w),
    Neg_share  = mean(w < 0),
    Pct_neg    = mean(w < 0) * 100,
    Sum_top10  = sum_top10,
    Sum_w      = sum(w),
    Sum_abs_w  = sum(abs(w)),
    Max_abs_w  = max(abs(w)),
    ESS_kish   = 1 / sw2,
    ESS_mod    = (sum(abs(w)))^2 / sw2
  )
}

# Format compact weight diagnostics.
weight_diag <- function(w, name) {
  stats <- weight_statistics(w)
  data.frame(
    Estimator = name,
    Sum_w     = round(stats[["Sum_w"]], 8),
    ESS       = round(stats[["ESS_kish"]], 0),
    ESS_mod   = round(stats[["ESS_mod"]], 0),
    Pct_neg   = round(stats[["Pct_neg"]], 1),
    Max_abs_w = round(stats[["Max_abs_w"]], 6),
    stringsAsFactors = FALSE
  )
}

# Calculate treatment-specific weight diagnostics.
outcome_weight_diagnostics <- function(w, D, name) {
  w <- as.numeric(w); D <- as.numeric(D)
  if (!length(w) || length(w) != length(D)) {
    stop("w and D must be nonempty vectors of equal length.")
  }
  if (any(!is.finite(w)) || any(!D %in% c(0, 1))) {
    return(data.frame(
      Estimator = name, Sum_w = NA_real_, Mass_T = NA_real_, Mass_C = NA_real_,
      ESS_mod_all = NA_real_, ESS_mod_T = NA_real_, ESS_mod_C = NA_real_,
      Wrong_mass_T = NA_real_, Wrong_mass_C = NA_real_, Max_abs_w = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  oriented <- ifelse(D == 1, w, -w)
  # Calculate modified ESS.
  ess_mod <- function(x) {
    den <- sum(x^2)
    if (length(x) == 0L || den == 0) return(NA_real_)
    sum(abs(x))^2 / den
  }
  # Calculate opposite-sign mass.
  wrong_mass <- function(x) {
    den <- sum(abs(x))
    if (length(x) == 0L || den == 0) return(NA_real_)
    sum(abs(x[x < 0])) / den
  }

  data.frame(
    Estimator = name,
    Sum_w = sum(w),
    Mass_T = sum(oriented[D == 1]),
    Mass_C = sum(oriented[D == 0]),
    ESS_mod_all = ess_mod(w),
    ESS_mod_T = ess_mod(w[D == 1]),
    ESS_mod_C = ess_mod(w[D == 0]),
    Wrong_mass_T = wrong_mass(oriented[D == 1]),
    Wrong_mass_C = wrong_mass(oriented[D == 0]),
    Max_abs_w = max(abs(w)),
    stringsAsFactors = FALSE
  )
}

# Format treatment-specific diagnostics.
outcome_weight_diag <- function(w, D, name) {
  x <- outcome_weight_diagnostics(w, D, name)
  # Format treated and control values.
  fmt_pair <- function(a, b, digits) {
    if (any(!is.finite(c(a, b)))) return(NA_character_)
    paste0(formatC(a, format = "f", digits = digits), " / ",
           formatC(b, format = "f", digits = digits))
  }
  data.frame(
    Estimator = x$Estimator,
    Sum_w = round(x$Sum_w, 8),
    `Mass T/C` = fmt_pair(x$Mass_T, x$Mass_C, 3),
    `ESS mod overall` = round(x$ESS_mod_all, 0),
    `ESS mod T/C` = fmt_pair(x$ESS_mod_T, x$ESS_mod_C, 0),
    `ESS mod/N T/C` = fmt_pair(x$ESS_mod_T / sum(D == 1),
                               x$ESS_mod_C / sum(D == 0), 3),
    `Opposite-sign mass T/C (%)` = fmt_pair(100 * x$Wrong_mass_T,
                                         100 * x$Wrong_mass_C, 1),
    `Max |w|` = round(x$Max_abs_w, 6),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# Calculate Kish ESS for positive weights.
positive_weight_ess <- function(w) {
  w <- as.numeric(w)
  if (!length(w) || any(!is.finite(w)) || any(w < 0) || sum(w^2) == 0) {
    return(NA_real_)
  }
  sum(w)^2 / sum(w^2)
}

# Build one weight-statistics block.
weight_stats_block <- function(w) {
  stats <- weight_statistics(w)
  c(
    n          = stats[["n"]],
    Min        = stats[["Min"]],
    Max        = stats[["Max"]],
    Neg_share  = round(stats[["Neg_share"]], 4),
    Sum_top10  = round(stats[["Sum_top10"]], 6),
    Sum_w      = round(stats[["Sum_w"]], 8),
    Sum_abs_w  = round(stats[["Sum_abs_w"]], 6),
    Max_abs_w  = round(stats[["Max_abs_w"]], 6),
    ESS_kish   = round(stats[["ESS_kish"]], 0),
    ESS_mod    = round(stats[["ESS_mod"]], 0)
  )
}

# Compare weight diagnostics by treatment status.
weight_bridge <- function(w, D, name) {
  w <- as.numeric(w); D <- as.numeric(D)
  blocks <- list(
    Overall = weight_stats_block(w),
    Treated = weight_stats_block(w[D == 1]),
    Control = weight_stats_block(-w[D == 0])
  )
  do.call(rbind, lapply(names(blocks), function(scope) {
    data.frame(Estimator = name, Scope = scope,
               as.list(blocks[[scope]]),
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
}

# Check an outcome-weight identity.
check_weight_identity <- function(w, Y, estimate) {
  isTRUE(all.equal(
    as.numeric(sum(w * Y)),
    as.numeric(estimate)
  ))
}

# Format estimates and standard errors.
fmt <- function(coef, se, digits = 3) {
  if (is.na(coef)) return("NA")
  if (is.na(se) || !is.finite(se) || se <= 0)
    return(sprintf(paste0("%.", digits, "f\n(NA)"), round(coef, digits)))
  pval  <- 2 * pnorm(-abs(coef / se))
  stars <- ifelse(pval < 0.01, "***", ifelse(pval < 0.05, "**",
                                             ifelse(pval < 0.10, "*", "")))
  sprintf(paste0("%.", digits, "f%s\n(%.", digits, "f)"),
          round(coef, digits), stars, round(se, digits))
}

# Create a numeric design matrix.
make_X <- function(formula, data) {
  mm  <- model.matrix(formula, data = data)
  m   <- matrix(as.numeric(mm), nrow = nrow(mm), ncol = ncol(mm))
  colnames(m) <- colnames(mm)
  m
}

# Extract a DML estimate.
get_estimate <- function(res, row_name) {
  if ("Estimate" %in% colnames(res)) {
    return(as.numeric(res[row_name, "Estimate"]))
  }
  as.numeric(res[row_name, 1])
}

# Check required outcome-weight rows.
check_omega_rows <- function(omega_obj, rows = c("PLR-IV", "Wald-AIPW")) {
  missing_rows <- setdiff(rows, rownames(omega_obj$omega))
  if (length(missing_rows) > 0) {
    stop("Missing expected outcome-weight rows: ",
         paste(missing_rows, collapse = ", "))
  }
}

# Extract a DoubleML coefficient.
get_dml_coef <- function(obj) {
  as.numeric(obj$coef)
}

# Extract a DoubleML standard error.
get_dml_se <- function(obj) {
  direct_se <- tryCatch(as.numeric(obj$se), error = function(e) numeric(0))
  if (length(direct_se) == 1L && is.finite(direct_se)) return(direct_se)

  s <- NULL
  invisible(capture.output(s <- obj$summary()))
  if ("std err"    %in% colnames(s)) return(as.numeric(s[1, "std err"]))
  if ("Std. Error" %in% colnames(s)) return(as.numeric(s[1, "Std. Error"]))
  NA_real_
}

# Construct all kappa weight vectors.
kappa_weights_bundle <- function(Z, D, X_kappa) {
  p_ml <- logit_mle(Z, X_kappa)
  p_cb <- get_cbps_p(Z, X_kappa)

  kw_ml <- kappa_outcome_weights(Z, D, p_ml)
  kw_cb <- kappa_outcome_weights(Z, D, p_cb)

  list(
    tau_cb_u   = kw_cb$w_u,
    tau_ml_u   = kw_ml$w_u,
    tau_ml_a10 = kw_ml$w_a10,
    tau_ml_a   = kw_ml$w_a,
    tau_ml_t   = kw_ml$w_a1,
    tau_ml_a0  = kw_ml$w_a0
  )
}

# Draw a covariate-balance Love plot.
make_love <- function(title_str, w_vec, D, X_bal, threshold = 0.1) {
  love.plot(
    D ~ X_bal,
    weights    = w_vec * (2 * D - 1),
    position   = "bottom",
    title      = title_str,
    thresholds = c(m = threshold),
    var.order  = "unadjusted",
    continuous = "std",
    binary     = "raw",
    s.d.denom  = "pooled",
    stats      = "mean.diffs",
    abs        = TRUE,
    line       = TRUE,
    stars      = "raw",
    colors     = viridis(2),
    shapes     = c("circle", "triangle")
  )
}

# Summarize a translation rerun.
translation_rerun_row <- function(name, tau_original, tau_shifted, k,
                                  tolerance = 1e-8, digits = 12) {
  difference <- as.numeric(tau_shifted) - as.numeric(tau_original)
  data.frame(
    Estimator        = name,
    Shift_k          = round(as.numeric(k), digits),
    tau_original     = round(as.numeric(tau_original), digits),
    tau_shifted      = round(as.numeric(tau_shifted), digits),
    rerun_difference = round(difference, digits),
    invariant        = is.finite(difference) && abs(difference) <= tolerance,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
}
