# ==============================================================================
# functions_all_v3.R
# ==============================================================================
# SINGLE shared function library for the thesis:
#   "Comparing Kappa Weighting and Causal Machine Learning Estimators
#    via Weight-Based Diagnostics"
#
# This is the ONLY function file you need. Source it once at the top of every
# replication script (Vietnam, Card, Angrist-Evans; kappa and DML):
#
#     source("functions_all_v3.R")
#
# It merges the former functions_kappa.R (kappa estimators, outcome weights,
# analytical SEs, 2SLS) with the DML / OutcomeWeights comparison helpers that
# were previously defined inline inside the DML notebooks. After sourcing this
# file, NO function needs to be defined inside any .Rmd.
#
# Variable conventions across all applications:
#   Y   : outcome variable
#   Z   : instrument (binary)
#   D   : treatment (binary)
#   p   : propensity score P(Z = 1 | X_i)
#   X   : covariate design matrix (including intercept column, for kappa)
#   X_df: covariate data frame (no intercept) for run_2sls()
#
# PACKAGES (load in the calling script, not here):
#   library(AER)            # ivreg() for run_2sls()
#   library(sandwich)       # vcovHC() for run_2sls()
#   library(lmtest)         # coeftest() for run_2sls()
#   library(OutcomeWeights) # dml_with_smoother(), get_outcome_weights()
#   library(DoubleML)       # DoubleMLIIVM etc. (extended DML analysis only)
#   library(mlr3verse)      # mlr3 learners (extended DML analysis only)
#   library(cobalt)         # love.plot()
#   library(viridis)        # colour palette
#   library(gridExtra)      # grid.arrange()
#
# ------------------------------------------------------------------------------
# SECTION MAP
# ------------------------------------------------------------------------------
# PART 1 — KAPPA ESTIMATORS (from former functions_kappa.R, unchanged)
#   §a  safe_logit()                — numerically stable logistic function
#   §b  prep_design_for_mest()      — standardise X for M-estimation
#   §c  logit_mle()                 — logit MLE propensity score (point est.)
#   §d  fit_logit_alpha()           — logit MLE for M-estimation sandwich
#   §e  fit_cbps_alpha()            — CBPS propensity score (Newton, line search)
#   §f  cbps() / get_cbps_p()       — CBPS convenience wrappers   [sig: (Z, X)]
#   §g  kappa_weights()             — kappa, kappa1, kappa0 (Abadie 2003)
#   §h  tau_u()                     — Uysal (2011) normalized estimator
#   §i  tau_a10()                   — Abadie-Cattaneo normalized estimator
#   §j  tau_unnorm()                — unnormalized kappa estimators (a/a1/a0)
#   §k  kappa_outcome_weights()     — omega_i weight vectors for all estimators
#   §l  num_jacobian()              — numerical Jacobian (central differences)
#   §m  matrix_inverse_safe()       — robust matrix inverse (fallback chain)
#   §n  sandwich_se_mest()          — sandwich SE from stacked moment function
#   §o  alpha_moment_matrix()       — propensity-score moment matrix (MLE/CBPS)
#   §p  kappa_analytic_se_one()     — M-estimation SE for one estimator
#   §q  safe_kappa_se()             — error-safe wrapper for kappa SE
#   §r  kappa_analytic_se_all()     — all six estimates + analytical SEs
#   §s  run_2sls()                  — 2SLS benchmark with HC1-robust SEs
#   §t  weight_diag()               — weight diagnostics (ESS, % neg, max|w|)
#   §u  check_weight_identity()     — algebraic check sum(w*Y) == tau_hat
#   §v  fmt()                       — coefficient formatter with stars
#                                     [UNUSED by Card notebook — see v3 note 4]
#
# PART 2 — DML / OUTCOMEWEIGHTS HELPERS (merged from the DML notebooks)
#   §A  make_X()                    — clean numeric design matrix for grf
#   §B  get_estimate()              — extract a point estimate from summary()
#   §C  check_omega_rows()          — assert expected rows exist in omega object
#   §D  check_doubleml_identity()   — omega'Y == coef check for DoubleML objects
#   §E  get_dml_coef()              — extract coef() from a DoubleML object
#   §E2 get_dml_se()                — extract SE from a DoubleML summary()
#                                     [PROMOTED from Card notebook in v3 — note 1]
#   §F  kappa_estimates()           — six kappa point estimates for one X spec
#                                     [UNUSED by Card notebook — see v3 note 4]
#   §G  kappa_weights_bundle()      — six kappa weight vectors for one X spec
#   §H  (removed in v3 — doubleml_translation_row() was dead — see note 2)
#   §I  make_love()                 — Love plot for one weight vector
#
# ------------------------------------------------------------------------------
# CONSISTENCY CHANGELOG (what changed when merging the notebook helpers)
# ------------------------------------------------------------------------------
# The DML notebooks previously redefined several functions inline, some of
# which silently OVERRODE the global kappa versions. Those local copies are
# removed. The single definitions below are the source of truth. Specifically:
#
#   • kappa_outcome_weights(), weight_diag(), check_weight_identity():
#       The notebooks redefined these locally with bodies identical to the
#       kappa versions. Now defined ONCE (§k, §t, §u). No behavioural change.
#
#   • get_cbps_p():  *** SIGNATURE FIX ***
#       The kappa version has signature get_cbps_p(Z, X) (§f).
#       The DML notebooks defined a one-argument get_cbps_p(X) that used Z as
#       a hidden global, and kappa_estimates()/kappa_weights_bundle() called
#       it as get_cbps_p(X_kappa). That one-argument form is REMOVED. The
#       two-argument §f version is the only one. kappa_estimates() (§F) and
#       kappa_weights_bundle() (§G) now call get_cbps_p(Z, X_kappa) correctly.
#
#   • kappa_estimates(), kappa_weights_bundle():  *** NO HIDDEN GLOBALS ***
#       Now take Z and D as explicit arguments (notebook used globals).
#
#   • doubleml_translation_row():  *** NO HIDDEN GLOBALS, DE-DUPLICATED ***
#       Was defined twice in the extended notebook and used Y_dol/Y_cnt/k as
#       globals. Now defined once (§H), takes Y_dollars/Y_cents/k explicitly.
#
#   • make_love():  *** NO HIDDEN GLOBALS, UNIFIED ***
#       Base notebook (make_love) and extended notebook (make_love_dml) both
#       hardcoded X_dml_cub and used D as a global. Unified into one make_love()
#       (§I) taking the balance matrix X_bal and treatment D as arguments.
#
#   • check_doubleml_identity() (§D) is INTENTIONALLY distinct from
#     check_weight_identity() (§u): the former takes a DoubleML omega OBJECT
#     (with an $omega matrix) and the object's coef; the latter takes a plain
#     weight vector and a scalar estimate. Both are kept.
#
# ------------------------------------------------------------------------------
# v2 -> v3 CLEANING PASS  (no econometrics changed; only consolidation + pruning)
# ------------------------------------------------------------------------------
# 1. PROMOTED get_dml_se() into the engine (§E2), directly after get_dml_coef().
#    It was the ONE helper still defined inline in card_presentation_3.Rmd
#    (its own comment admitted "the only helper not in functions_all.R").
#    Body copied verbatim. After this, the inline definition can be deleted
#    from the .Rmd and the header's "NO function defined inside any .Rmd"
#    claim is finally true. No behavioural change.
#
# 2. REMOVED doubleml_translation_row() (was §H). It was DEAD CODE: not called
#    anywhere in card_presentation_3.Rmd (it survived only in a stale comment),
#    and it is a strict subset of ti_unified_row() (§K), which computes the same
#    cents-vs-dollars Method-B shift plus the Method-A rerun in one row.
#    >>> If your Vietnam notebook still calls doubleml_translation_row(), tell me
#        and I will restore it in one line. Nothing in the Card pipeline uses it.
#
# 3. NAMING INCONSISTENCY — FLAGGED, NOT FIXED (fixing it now would break the
#    Card notebook). The Tan/Frolich unnormalized estimator (kappa1 denominator)
#    is keyed differently across functions:
#         kappa_analytic_se_all()  (§r)  ->  "tau_ml_t"
#         kappa_estimates()        (§F)  ->  "tau_ml_a1"
#         kappa_weights_bundle()   (§G)  ->  "tau_ml_a1"
#    The Card notebook relies on BOTH spellings in different chunks, so renaming
#    here would break it. This should be unified (one canonical key) when we
#    next touch the notebooks — ask me and we do it across engine + notebooks
#    together.
#
# 4. UNUSED-BY-CARD (kept, not removed — flagged inline with an [UNUSED] tag):
#    fmt() (§v), kappa_estimates() (§F), translation_check_weights() and
#    translation_check_rerun() (§J). None is called by card_presentation_3.Rmd.
#    They are retained because (a) fmt() is plausibly used by the Vietnam
#    notebook, (b) the two standalone translation_check_* functions are the
#    natural building blocks for your planned SEPARATE translation-invariance
#    notebook, and (c) kappa_estimates() is a lighter no-SE alternative to
#    kappa_analytic_se_all()$estimates. Say the word and any of these can go.
#
# 5. NOT TOUCHED: the econometric bodies (propensity scores, kappa formulas,
#    sandwich SEs, outcome-weight constructors) are byte-for-byte unchanged.
#    This pass did not re-verify the math; it only consolidated and pruned.
# ==============================================================================


# ==============================================================================
# ==============================================================================
# PART 1 — KAPPA ESTIMATORS, OUTCOME WEIGHTS, ANALYTICAL SEs, 2SLS
# ==============================================================================
# ==============================================================================
# These functions are reproduced verbatim from the former functions_kappa.R.
# They are unchanged; only the surrounding file and the DML helpers (Part 2)
# are new.


# ==============================================================================
# §a  NUMERICALLY STABLE LOGISTIC FUNCTION
# ==============================================================================
# Clips eta to [-35, 35] to prevent exp() overflow or underflow.
# Used in all propensity score computations.

safe_logit <- function(eta) {
  eta <- pmin(pmax(eta, -35), 35)
  1 / (1 + exp(-eta))
}


# ==============================================================================
# §b  DESIGN MATRIX STANDARDISATION
# ==============================================================================
# Centres and scales non-intercept columns for better numerical conditioning
# of the Jacobian in M-estimation. The intercept column (all 1s) is left
# unchanged. Fitted propensity scores are identical before and after
# standardisation (pure reparameterisation).

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


# ==============================================================================
# §c  LOGIT MLE PROPENSITY SCORE  [point estimates only]
# ==============================================================================
# Standard glm() logit. The "-1" suppresses a duplicate intercept when X
# already contains a column of 1s. Returns fitted P(Z = 1 | X_i).
# Used whenever point estimates only are needed (not M-estimation SEs).

logit_mle <- function(Z, X) {
  df  <- data.frame(Z = Z, X)
  fit <- glm(Z ~ . - 1, data = df, family = binomial(link = "logit"))
  fitted.values(fit)
}


# ==============================================================================
# §d  LOGIT MLE FOR M-ESTIMATION  [returns alpha + standardised X]
# ==============================================================================
# Used internally by kappa_analytic_se_one(). Standardises X via
# prep_design_for_mest() for numerical stability of the Jacobian.
# Returns alpha (on standardised scale), fitted p, and X_used.

fit_logit_alpha <- function(Z, X) {
  Z     <- as.numeric(Z)
  X     <- prep_design_for_mest(X)
  fit   <- glm.fit(x = X, y = Z, family = binomial(link = "logit"))
  alpha <- as.numeric(coef(fit));  alpha[is.na(alpha)] <- 0
  p     <- as.vector(pmin(pmax(safe_logit(X %*% alpha), 1e-8), 1 - 1e-8))
  list(alpha = alpha, p = p, X_used = X)
}


# ==============================================================================
# §e  CBPS PROPENSITY SCORE  [Newton with backtracking line search]
# ==============================================================================
# Solves the covariate-balancing moment condition
#
#   E[ (Z - p(X)) / {p(X)(1-p(X))} * X ] = 0
#
# via Newton steps with backtracking. Initialised at the logit MLE.
# Returns the best iterate (lowest max-moment norm) if tol is not reached.
# Used exclusively for tau_cb_u (CBPS + Uysal normalisation).

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

  moment_fn <- function(b) {
    p <- as.vector(pmin(pmax(safe_logit(X %*% b), 1e-8), 1 - 1e-8))
    colMeans(as.vector((Z - p) / (p * (1 - p))) * X)
  }

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


# ==============================================================================
# §f  CBPS CONVENIENCE WRAPPERS
# ==============================================================================
# cbps()      : calls fit_cbps_alpha(), returns the full list.
# get_cbps_p(): calls cbps() and extracts only the p vector.
#   Used wherever a plain propensity score vector is needed.
#
# CHANGELOG: get_cbps_p() was absent in vietnam_14_05.R (CBPS p was extracted
# inline). Promoted here from card_21_05.R — no behavioural change.

cbps <- function(Z, X, tol = 1e-9, max_iter = 5000, verbose = FALSE) {
  fit_cbps_alpha(Z, X, tol = tol, max_iter = max_iter)
}

get_cbps_p <- function(Z, X) {
  out <- cbps(Z, X)
  if (is.list(out) && !is.null(out$p)) return(as.vector(out$p))
  as.vector(out)
}


# ==============================================================================
# §g  KAPPA WEIGHTS  (Abadie 2003, Lemma 2.1)
# ==============================================================================
# The three kappa weights identify complier moments:
#
#   kappa  = 1 - D(1-Z)/(1-p) - (1-D)Z/p
#   kappa1 = D(Z - p) / [p(1-p)]
#   kappa0 = (1-D)((1-Z)-(1-p)) / [p(1-p)]
#
# In population: E[kappa] = E[kappa1] = E[kappa0] = P(complier).

kappa_weights <- function(Z, D, p) {
  list(
    kappa  = 1 - D * (1 - Z) / (1 - p) - (1 - D) * Z / p,
    kappa1 = D * (Z - p) / (p * (1 - p)),
    kappa0 = (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
  )
}


# ==============================================================================
# §h  tau_u — UYSAL (2011) NORMALIZED ESTIMATOR  [translation invariant]
# ==============================================================================
# Computes separately normalised IPW means for Z=1 and Z=0, then takes
# the ratio of differences:
#
#   tau_u = [mu_Y1 - mu_Y0] / [mu_D1 - mu_D0]
#
# Translation invariant: adding a constant c to Y shifts both mu_Y1 and
# mu_Y0 by c, so the numerator (and hence tau_u) is unchanged.

tau_u <- function(Y, Z, D, p) {
  s1 <- sum(Z / p)
  s0 <- sum((1 - Z) / (1 - p))
  numerator   <- sum(Y * Z / p) / s1 - sum(Y * (1 - Z) / (1 - p)) / s0
  denominator <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  numerator / denominator
}


# ==============================================================================
# §i  tau_a10 — ABADIE-CATTANEO NORMALIZED ESTIMATOR  [translation invariant]
# ==============================================================================
# Separately normalises kappa1 and kappa0 weighted outcome means:
#
#   tau_a10 = sum(kappa1 * Y) / sum(kappa1)
#           - sum(kappa0 * Y) / sum(kappa0)
#
# Translation invariant because the constant c cancels in each ratio.

tau_a10 <- function(Y, Z, D, p) {
  kw <- kappa_weights(Z, D, p)
  sum(kw$kappa1 * Y) / sum(kw$kappa1) - sum(kw$kappa0 * Y) / sum(kw$kappa0)
}


# ==============================================================================
# §j  tau_unnorm — UNNORMALIZED KAPPA ESTIMATORS  [NOT translation invariant]
# ==============================================================================
# Common numerator: Delta = mean[ Y * (Z - p) / {p(1-p)} ]
# Three denominator choices (SUW notation):
#   "a"  : Gamma = mean(kappa)    — Abadie (2003) original
#   "a1" : Gamma = mean(kappa1)   — Tan (2006) / Frölich (2007); tau_t in SUW
#   "a0" : Gamma = mean(kappa0)
#
# NOT translation invariant: adding c to Y shifts Delta but not Gamma.

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


# ==============================================================================
# §k  OUTCOME-WEIGHT CONSTRUCTORS  (omega_i representation)
# ==============================================================================
# Every kappa estimator can be written as tau_hat = sum_i omega_i * Y_i.
# Key properties:
#   sum(w) == 0  ⟺  translation invariant  (holds for w_u and w_a10)
#   ESS = 1 / sum(w^2)
#   % negative weights = diagnostic for instability
#
# CHANGELOG: absent in vietnam_14_05.R; added in card_21_05.R; copied verbatim.

kappa_outcome_weights <- function(Z, D, p) {
  n  <- length(Z)
  kw <- kappa_weights(Z, D, p)

  # tau_u weights
  s1  <- sum(Z / p)
  s0  <- sum((1 - Z) / (1 - p))
  dD  <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  w_u <- (Z / p / s1 - (1 - Z) / (1 - p) / s0) / dD

  # tau_a10 weights
  w_a10 <- kw$kappa1 / sum(kw$kappa1) - kw$kappa0 / sum(kw$kappa0)

  # common numerator weight for unnormalized estimators
  num_w <- (Z - p) / (p * (1 - p)) / n

  list(
    w_u   = as.vector(w_u),
    w_a10 = as.vector(w_a10),
    w_a   = as.vector(num_w / mean(kw$kappa)),
    w_a1  = as.vector(num_w / mean(kw$kappa1)),
    w_a0  = as.vector(num_w / mean(kw$kappa0))
  )
}


# ==============================================================================
# §l  NUMERICAL JACOBIAN  (central differences)
# ==============================================================================
# Used inside sandwich_se_mest() to avoid hand-coding the Jacobian of every
# stacked moment system. Adaptive step size: h = eps * (|theta_j| + 1).
#
# CHANGELOG: vietnam_14_05.R used variable names th_plus/th_minus.
#            card_21_05.R used tp/tm. Global version uses tp/tm (Card
#            convention). Behaviour is identical.

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


# ==============================================================================
# §m  SAFE MATRIX INVERSE  (fallback chain)
# ==============================================================================
# Tries solve() → qr.solve() → ridge-regularised solve() → SVD pseudo-inverse.
# Returns the first finite result found.

matrix_inverse_safe <- function(A, tol = 1e-10) {
  inv <- tryCatch(solve(A), error = function(e) NULL)
  if (!is.null(inv) && all(is.finite(inv))) return(inv)

  inv <- tryCatch(qr.solve(A), error = function(e) NULL)
  if (!is.null(inv) && all(is.finite(inv))) return(inv)

  for (ridge in c(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)) {
    inv <- tryCatch(solve(A + ridge * diag(ncol(A))), error = function(e) NULL)
    if (!is.null(inv) && all(is.finite(inv))) return(inv)
  }

  # Moore-Penrose via SVD
  sv    <- svd(A)
  d     <- sv$d
  d_inv <- ifelse(d > tol * max(d), 1 / d, 0)
  sv$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(sv$u)
}


# ==============================================================================
# §n  SANDWICH SE FROM STACKED MOMENT FUNCTION  (M-estimation)
# ==============================================================================
# Computes the SE of theta[tau_index] from the sandwich formula
#
#   Avar(sqrt(n) * theta_hat) = A^{-1} V (A^{-1})'
#
# where:
#   A = E[d/dtheta' psi(O_i, theta)]  — Jacobian of mean moment (numerical)
#   V = Var(psi(O_i, theta_0))        — variance of moment contributions
#
# The SE for tau_index is sqrt(vcov[tau_index, tau_index]).

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


# ==============================================================================
# §o  PROPENSITY-SCORE MOMENT MATRIX
# ==============================================================================
# Returns the n x k matrix of per-observation propensity-score moment
# contributions psi_alpha_i, stacked with the LATE estimator moments.
#
#   logit MLE: psi_alpha_i = (Z_i - p_i) * X_i   [score of logit log-lik]
#   CBPS:      psi_alpha_i = (Z_i - p_i) / {p_i(1-p_i)} * X_i
#
# CHANGELOG: absent in vietnam_14_05.R (was inline inside kappa_analytic_se_one);
#            extracted into its own function in card_21_05.R; copied verbatim.

alpha_moment_matrix <- function(Z, p, X_used, method) {
  p <- as.vector(p);  Z <- as.vector(Z)
  if (method == "ml") return(as.vector(Z - p) * X_used)
  if (method == "cb") return(as.vector((Z - p) / (p * (1 - p))) * X_used)
  stop("method must be 'ml' or 'cb'")
}


# ==============================================================================
# §p  ANALYTICAL M-ESTIMATION SE FOR ONE KAPPA ESTIMATOR
# ==============================================================================
# Stacks the propensity-score moments with the LATE estimator moments and calls
# sandwich_se_mest(). The tau_index is the last component of theta_hat.
#
# Supported estimators:
#   "u"   — tau_u  (Uysal normalized)
#     theta = (alpha, mu1, mu0, m1, m0, tau), dim = k+5
#   "a10" — tau_a10
#     theta = (alpha, Delta1, Gamma1, Delta0, Gamma0, tau), dim = k+5
#   "a"   — tau_a  (Abadie unnormalized, kappa denominator)
#   "a1"  — tau_a1 = tau_t (Tan/Frölich)
#   "a0"  — tau_a0
#     theta = (alpha, Delta, Gamma, tau), dim = k+3

kappa_analytic_se_one <- function(Y, Z, D, X, estimator, method = "ml") {
  Y <- as.numeric(Y);  Z <- as.numeric(Z)
  D <- as.numeric(D);  X <- as.matrix(X)
  n <- length(Y)

  fit       <- if (method == "ml") fit_logit_alpha(Z, X) else fit_cbps_alpha(Z, X)
  alpha_hat <- fit$alpha
  X_used    <- fit$X_used
  k         <- length(alpha_hat)

  am <- function(p) alpha_moment_matrix(Z, p, X_used, method)

  # ---------------------------------------------------------------------------
  # tau_u: theta = (alpha [k], mu1, mu0, m1, m0, tau) — length k+5
  # ---------------------------------------------------------------------------
  if (estimator == "u") {
    p   <- fit$p
    s1  <- sum(Z / p);  s0 <- sum((1 - Z) / (1 - p))
    mu1 <- sum(Z * Y / p) / s1;            mu0 <- sum((1 - Z) * Y / (1 - p)) / s0
    m1  <- sum(Z * D / p) / s1;            m0  <- sum((1 - Z) * D / (1 - p)) / s0
    tau <- (mu1 - mu0) / (m1 - m0)

    theta_hat <- c(alpha_hat, mu1, mu0, m1, m0, tau)

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

  # ---------------------------------------------------------------------------
  # tau_a, tau_a1, tau_a0: theta = (alpha [k], Delta, Gamma, tau) — length k+3
  # ---------------------------------------------------------------------------
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

  # ---------------------------------------------------------------------------
  # tau_a10: theta = (alpha [k], Delta1, Gamma1, Delta0, Gamma0, tau) — k+5
  # ---------------------------------------------------------------------------
  if (estimator == "a10") {
    p   <- fit$p
    kw  <- kappa_weights(Z, D, p)
    Delta1 <- mean(kw$kappa1 * Y);  Gamma1 <- mean(kw$kappa1)
    Delta0 <- mean(kw$kappa0 * Y);  Gamma0 <- mean(kw$kappa0)
    tau    <- Delta1/Gamma1 - Delta0/Gamma0
    theta_hat <- c(alpha_hat, Delta1, Gamma1, Delta0, Gamma0, tau)

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


# ==============================================================================
# §q  ERROR-SAFE WRAPPER FOR KAPPA SE
# ==============================================================================
# Returns NA_real_ (with a warning) if the sandwich system is singular,
# rather than stopping the script mid-table.

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


# ==============================================================================
# §r  FULL KAPPA TABLE: ESTIMATES + ANALYTICAL SEs
# ==============================================================================
# Convenience wrapper that computes all six kappa point estimates and their
# corresponding M-estimation standard errors in one call.
#
# Returns:
#   list(estimates = named numeric(6), se = named numeric(6))
#
# Estimator names (matching SUW Table notation):
#   tau_cb_u   : CBPS + Uysal normalisation      [Panel B, row 1]
#   tau_ml_u   : MLE  + Uysal normalisation      [Panel B, row 2]
#   tau_ml_a10 : MLE  + Abadie-Cattaneo          [Panel B, row 3]
#   tau_ml_a   : MLE  + unnorm, kappa denom      [Panel C, row 4]
#   tau_ml_t   : MLE  + unnorm, kappa1 denom     [Panel C, row 5]
#   tau_ml_a0  : MLE  + unnorm, kappa0 denom     [Panel C, row 6]

kappa_analytic_se_all <- function(Y, Z, D, X_mat) {
  Y <- as.numeric(Y); Z <- as.numeric(Z)
  D <- as.numeric(D); X_mat <- as.matrix(X_mat)

  p_ml <- logit_mle(Z, X_mat)
  p_cb <- get_cbps_p(Z, X_mat)

  estimates <- c(
    tau_cb_u   = tau_u(Y, Z, D, p_cb),
    tau_ml_u   = tau_u(Y, Z, D, p_ml),
    tau_ml_a10 = tau_a10(Y, Z, D, p_ml),
    tau_ml_a   = tau_unnorm(Y, Z, D, p_ml, "a"),
    tau_ml_t   = tau_unnorm(Y, Z, D, p_ml, "a1"),   # tau_t = tau_a,1 in SUW
                                                     # FLAG: keyed "tau_ml_a1" in
                                                     # §F/§G — unify later (note 3)
    tau_ml_a0  = tau_unnorm(Y, Z, D, p_ml, "a0")
  )

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


# ==============================================================================
# §s  2SLS BENCHMARK WITH HC1-ROBUST SEs
# ==============================================================================
# Standard ivreg() from the AER package with sandwich HC1 robust SEs.
# HC1 = (n/(n-k)) * HC0 — matches Stata's ivreg2 vce(robust) default.
# Formula: Y ~ D + X_controls | Z + X_controls
#
# Arguments:
#   Y          : outcome vector
#   D          : treatment vector
#   Z          : instrument vector
#   X_df       : covariate data.frame (NO intercept column)
#   endog_name : name for D in the formula (default "D")

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


# ==============================================================================
# §t  WEIGHT DIAGNOSTICS
# ==============================================================================
# Key diagnostics for a weight vector omega:
#   Sum_w     : sum(omega_i) — should be ~0 for translation-invariant estimators
#   ESS       : Kish (1965) effective sample size = 1 / sum(omega_i^2).
#               *** Only interpretable for NONNEGATIVE, sum-to-one weights. ***
#               For outcome weights (which carry negatives and do NOT sum to 1)
#               it can fall below 1 or exceed n, so it is kept for continuity
#               but should NOT be the reported headline — use ESS_mod instead.
#   ESS_mod   : modified ESS = (sum|omega_i|)^2 / sum(omega_i^2)
#               Chattopadhyay & Zubizarreta (2023, Biometrika 110(3), §5.4,
#               doi:10.1093/biomet/asac058). Scale-invariant, guaranteed in
#               [1, n], and reduces EXACTLY to the Kish ESS when all weights are
#               nonnegative. This is the effective-sample-size number to report
#               and to compare across kappa / DML estimators.
#   Pct_neg   : % observations with omega_i < 0
#   Max_abs_w : max |omega_i| — outlier detection
#
# v4: added ESS_mod (the modified ESS). This widens the returned data.frame
#     from 5 to 6 columns; extend any align=c("l","r","r","r","r") vector in the
#     notebooks to 6 entries (cosmetic only — knitting is unaffected).
# CHANGELOG: absent in vietnam_14_05.R; added in card_21_05.R.

weight_diag <- function(w, name) {
  sw2 <- sum(w^2)
  data.frame(
    Estimator = name,
    Sum_w     = round(sum(w), 8),
    ESS       = round(1 / sw2, 0),                    # Kish — see caveat above
    ESS_mod   = round((sum(abs(w)))^2 / sw2, 0),      # Zubizarreta modified ESS
    Pct_neg   = round(mean(w < 0) * 100, 1),
    Max_abs_w = round(max(abs(w)), 6),
    stringsAsFactors = FALSE
  )
}



# ==============================================================================
# §t2  WEIGHT-STATISTICS BRIDGE  (RQ1: my diagnostics  <->  Knaus's summary())
# ==============================================================================
# Computes, for ONE outcome-weight vector, the union of the two diagnostic
# vocabularies used in the thesis, reported Overall / Treated / Control:
#
#   (a) "Knaus block" — the six statistics his summary.get_outcome_weights()
#       prints, reimplemented here in plain R so they are available for the
#       KAPPA estimators too (his package only produces them for DoubleML
#       objects). Definitions match his C++ summary_weights_rcpp() exactly:
#         Min, Max, Pct_neg, Sum_top10 (sum of the largest 10% of weights),
#         Sum_w, Sum_abs_w.
#       Group convention also matches his summary(): Treated uses omega[D==1]
#       as-is; Control uses the SIGN-FLIPPED -omega[D==0], so control weights
#       are summarised on the same positive scale he uses.
#
#   (b) "ESS block" — my effective-sample-size measures:
#         ESS_kish = 1 / sum(w^2)                 (Kish 1965; caveat in §t)
#         ESS_mod  = (sum|w|)^2 / sum(w^2)         (Zubizarreta 2023, §5.4)
#       ESS_mod is sign-flip invariant, so Treated/Control values are the same
#       whether or not the control sign flip is applied.
#
# This single object is the computational core of RQ1: it puts "what he defines"
# and "what I define" in one table and lets you show they describe the same
# weights. It does NOT replace weight_diag() (§t) — that stays as the compact
# per-estimator row used in the existing notebook tables.
#
# Arguments:
#   w    : signed outcome-weight vector (length n)
#   D    : binary treatment indicator (length n) — defines Treated/Control split
#   name : estimator label
#
# Returns: data.frame, 3 rows (Scope = Overall / Treated / Control), columns
#   Estimator, Scope, n, Min, Max, Pct_neg, Sum_top10, Sum_w, Sum_abs_w,
#   ESS_kish, ESS_mod.

weight_stats_block <- function(w) {
  w   <- as.numeric(w)
  n   <- length(w)
  sw2 <- sum(w^2)
  ord <- sort(w)
  p90 <- floor(0.9 * n)                              # matches his C++ index
  sum_top10 <- if (p90 < n) sum(ord[(p90 + 1):n]) else sum(ord)
  c(
    n          = n,
    Min        = min(w),
    Max        = max(w),
    Pct_neg    = round(mean(w < 0) * 100, 1),
    Sum_top10  = round(sum_top10, 6),
    Sum_w      = round(sum(w), 8),
    Sum_abs_w  = round(sum(abs(w)), 6),
    ESS_kish   = round(1 / sw2, 0),
    ESS_mod    = round((sum(abs(w)))^2 / sw2, 0)
  )
}

weight_bridge <- function(w, D, name) {
  w <- as.numeric(w); D <- as.numeric(D)
  blocks <- list(
    Overall = weight_stats_block(w),
    Treated = weight_stats_block(w[D == 1]),
    Control = weight_stats_block(-w[D == 0])         # Knaus sign-flip convention
  )
  do.call(rbind, lapply(names(blocks), function(scope) {
    data.frame(Estimator = name, Scope = scope,
               as.list(blocks[[scope]]),
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
}



# ==============================================================================
# §u  ALGEBRAIC IDENTITY CHECK
# ==============================================================================
# Verifies: sum(w * Y) == tau_hat  (should hold to machine precision).
# Returns TRUE/FALSE. Use after computing outcome weights to confirm that
# the weight vector correctly represents the estimator.
#
# CHANGELOG: absent in vietnam_14_05.R; added in card_21_05.R; copied verbatim.

check_weight_identity <- function(w, Y, estimate, tol = 1e-8) {
  isTRUE(all.equal(sum(w * Y), estimate, tolerance = tol))
}


# ==============================================================================
# §v  OUTPUT FORMATTER     [UNUSED by card_presentation_3.Rmd — kept for Vietnam]
# ==============================================================================
# Formats a coefficient with significance stars and standard error in
# parentheses. Stars: *** p<0.01, ** p<0.05, * p<0.10 (two-sided z-test).
# Returns a character string for cat() or table output.
# NOTE: the Card notebook uses its own inline cell formatters (cell_fmt / cell /
# ext_cell) instead of fmt(). Retained because the Vietnam notebook likely calls
# it; remove if you confirm nothing uses it.

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


# ==============================================================================
# ==============================================================================
# PART 2 — DML / OUTCOMEWEIGHTS HELPERS
# ==============================================================================
# ==============================================================================
# All helpers below take every input as an explicit argument (no hidden
# globals). They depend on the PART 1 kappa functions above (logit_mle,
# get_cbps_p, tau_u, tau_a10, tau_unnorm, kappa_outcome_weights) and on the
# OutcomeWeights / DoubleML / cobalt packages loaded by the calling script.


# ==============================================================================
# §A  CLEAN NUMERIC DESIGN MATRIX FOR grf
# ==============================================================================
# Strips all attributes from model.matrix() output to produce a plain numeric
# matrix that passes grf's validate_X() checks. grf in R 4.5.x rejects
# single-column matrices, so the caller must ensure ncol >= 2.
# CONSISTENCY: identical to make_X() in the base Vietnam DML notebook.

make_X <- function(formula, data) {
  mm  <- model.matrix(formula, data = data)
  m   <- matrix(as.numeric(mm), nrow = nrow(mm), ncol = ncol(mm))
  colnames(m) <- colnames(mm)
  m
}


# ==============================================================================
# §B  EXTRACT A POINT ESTIMATE FROM summary(dml_object)
# ==============================================================================
# Robust to whether the summary names its first column "Estimate" or uses
# positional indexing. Used to pull "PLR-IV" / "Wald-AIPW" rows.
# CONSISTENCY: identical to get_estimate() in the base Vietnam DML notebook.

get_estimate <- function(res, row_name) {
  if ("Estimate" %in% colnames(res)) {
    return(as.numeric(res[row_name, "Estimate"]))
  }
  as.numeric(res[row_name, 1])
}


# ==============================================================================
# §C  ASSERT EXPECTED ROWS EXIST IN AN omega OBJECT
# ==============================================================================
# Stops immediately if get_outcome_weights() did not return the expected
# estimator rows (guards against package-version changes in row order/naming).
# CONSISTENCY: identical to check_omega_rows() in the base Vietnam DML notebook.

check_omega_rows <- function(omega_obj, rows = c("PLR-IV", "Wald-AIPW")) {
  missing_rows <- setdiff(rows, rownames(omega_obj$omega))
  if (length(missing_rows) > 0) {
    stop("Missing expected outcome-weight rows: ",
         paste(missing_rows, collapse = ", "))
  }
}


# ==============================================================================
# §D  ALGEBRAIC IDENTITY CHECK FOR A DoubleML OBJECT
# ==============================================================================
# DoubleML analogue of check_weight_identity() (§u). Verifies omega'Y == coef
# to machine precision for a DoubleML-derived omega object. Returns TRUE/FALSE.
# This is the Pathway-2 diagnostic: TRUE for lm/ranger, typically FALSE (~1e-6)
# for XGBoost, reflecting the affine-smoother caveat in Knaus (2024).
# CONSISTENCY: identical to check_doubleml_identity() in the extended notebook.
# NOTE: distinct from §u — that takes a plain weight vector + scalar estimate;
# this takes a DoubleML omega OBJECT (with $omega matrix) + the object's coef.

check_doubleml_identity <- function(omega_obj, Y, coef, tol = 1e-8) {
  isTRUE(all.equal(
    as.numeric(omega_obj$omega %*% Y),
    as.numeric(coef),
    tolerance = tol
  ))
}


# ==============================================================================
# §E  EXTRACT coef() FROM A DoubleML OBJECT
# ==============================================================================
# Thin accessor so the comparison tables read uniformly.
# CONSISTENCY: identical to get_dml_coef() in the extended notebook.

get_dml_coef <- function(obj) {
  as.numeric(obj$coef)
}


# ==============================================================================
# §E2  EXTRACT SE FROM A DoubleML summary()
# ==============================================================================
# SE analogue of get_dml_coef() (§E). Pulls the standard error from a fitted
# DoubleML object's summary(), tolerating either column name used across
# DoubleML versions ("std err" or "Std. Error"). Returns NA_real_ if neither
# column is present.
#
# v3: PROMOTED verbatim from card_presentation_3.Rmd, where it was the only
# helper still defined inline. The comparison/learner tables (Linear+Logit,
# ranger, XGBoost) call this for every cell.

get_dml_se <- function(obj) {
  s <- obj$summary()
  if ("std err"    %in% colnames(s)) return(as.numeric(s[1, "std err"]))
  if ("Std. Error" %in% colnames(s)) return(as.numeric(s[1, "Std. Error"]))
  NA_real_
}


# ==============================================================================
# §F  SIX KAPPA POINT ESTIMATES FOR ONE COVARIATE SPECIFICATION
#                          [UNUSED by card_presentation_3.Rmd — see note below]
# ==============================================================================
# All six kappa estimates (CBPS+Uysal, MLE+Uysal, Abadie-Cattaneo, and the
# three unnormalised variants) for a given kappa design matrix X_kappa.
# tau_ml_a1 is tau_t in SUW notation. Returns a named numeric(6).
#
# UNUSED: the Card notebook reads its six point estimates from
# kappa_analytic_se_all()$estimates (§r) instead. This function is the lighter
# no-SE alternative; kept for convenience.
#
# *** NAMING INCONSISTENCY (flagged, not fixed — see v3 header note 3) ***
# This function keys the Tan/Frolich estimator as "tau_ml_a1", but
# kappa_analytic_se_all() (§r) keys the SAME estimator as "tau_ml_t". Do not
# rename here in isolation — the Card notebook depends on both spellings. Unify
# across engine + notebooks in one coordinated pass.
#
# *** CHANGED vs notebook ***  Z and D are explicit arguments (were globals),
# and get_cbps_p() is called with the correct two-argument signature
# get_cbps_p(Z, X_kappa). The notebook's one-argument get_cbps_p(X_kappa) only
# worked against a local one-arg redefinition that is now removed.

kappa_estimates <- function(Y, Z, D, X_kappa) {
  p_ml <- logit_mle(Z, X_kappa)
  p_cb <- get_cbps_p(Z, X_kappa)          # CHANGED: was get_cbps_p(X_kappa)

  c(
    tau_cb_u   = tau_u(Y, Z, D, p_cb),
    tau_ml_u   = tau_u(Y, Z, D, p_ml),
    tau_ml_a10 = tau_a10(Y, Z, D, p_ml),
    tau_ml_a   = tau_unnorm(Y, Z, D, p_ml, "a"),
    tau_ml_a1  = tau_unnorm(Y, Z, D, p_ml, "a1"),   # FLAG: keyed "tau_ml_t" in §r
    tau_ml_a0  = tau_unnorm(Y, Z, D, p_ml, "a0")
  )
}


# ==============================================================================
# §G  SIX KAPPA OUTCOME-WEIGHT VECTORS FOR ONE COVARIATE SPECIFICATION
# ==============================================================================
# Per-observation outcome weights for all six kappa estimators, via
# kappa_outcome_weights() (§k). Returns a named list of six weight vectors,
# matching the keys used by weight_diag() (§t) in the diagnostics tables.
#
# *** CHANGED vs notebook ***  Same fix as §F: Z and D explicit; get_cbps_p()
# called as get_cbps_p(Z, X_kappa).

kappa_weights_bundle <- function(Z, D, X_kappa) {
  p_ml <- logit_mle(Z, X_kappa)
  p_cb <- get_cbps_p(Z, X_kappa)          # CHANGED: was get_cbps_p(X_kappa)

  kw_ml <- kappa_outcome_weights(Z, D, p_ml)
  kw_cb <- kappa_outcome_weights(Z, D, p_cb)

  list(
    tau_cb_u   = kw_cb$w_u,
    tau_ml_u   = kw_ml$w_u,
    tau_ml_a10 = kw_ml$w_a10,
    tau_ml_a   = kw_ml$w_a,
    tau_ml_a1  = kw_ml$w_a1,
    tau_ml_a0  = kw_ml$w_a0
  )
}


# ==============================================================================
# §H  (REMOVED in v3)
# ==============================================================================
# doubleml_translation_row() lived here. It was never called by
# card_presentation_3.Rmd and is fully superseded by ti_unified_row() (§K),
# which reports the same cents-vs-dollars Method-B shift alongside the Method-A
# rerun. See the v3 cleaning note 2 in the header. To restore, ask.


# ==============================================================================
# §I  LOVE PLOT FOR ONE WEIGHT VECTOR
# ==============================================================================
# cobalt Love plot of absolute standardised mean differences, with the
# (2*D - 1) sign flip mapping signed outcome weights onto cobalt's
# treated-vs-control convention (Knaus 2024, 401k notebook).
#
# *** CHANGED vs notebook ***  Base (make_love) and extended (make_love_dml)
# both hardcoded X_dml_cub and used D as a global. Unified here: balance matrix
# X_bal and treatment D are explicit arguments, so it serves any specification
# (cubic, saturated, poly) and any application.
#
# Arguments:
#   title_str : panel title
#   w_vec     : signed outcome weight vector
#   D         : treatment indicator (for the 2D-1 flip and the formula LHS)
#   X_bal     : covariate matrix to assess balance on
#   threshold : ASMD reference line (default 0.1)

make_love <- function(title_str, w_vec, D, X_bal, threshold = 0.1) {
  love.plot(
    D ~ X_bal,
    weights    = w_vec * (2 * D - 1),
    position   = "bottom",
    title      = title_str,
    thresholds = c(m = threshold),
    var.order  = "unadjusted",
    binary     = "std",
    stats      = c("mean.diffs", "ks.statistics"), # stats = "mean.diffs"
    abs        = TRUE,
    line       = TRUE,
    colors     = viridis(2),
    shapes     = c("circle", "triangle")
  )
}


# ==============================================================================
# §J  TRANSLATION INVARIANCE — TWO METHODS, HIGH PRECISION
#     [translation_check_weights / translation_check_rerun are UNUSED by the
#      Card notebook — it uses the combined ti_unified_row() (§K) instead. These
#      standalone Method-A / Method-B builders are kept on purpose: they are the
#      natural pieces for your planned SEPARATE translation-invariance notebook,
#      where refitting in isolation is the whole point.]
# ==============================================================================
# Both methods test the same property: tau_hat(Y + c) = tau_hat(Y).
# Adding a constant c changes a weighted-sum estimator by exactly c * sum(w),
# so the estimator is translation invariant  <=>  sum_i w_i = 0.
#
# In a LOG-OUTCOME application, switching the wage units from dollars to cents
# multiplies the raw wage by 100, which ADDS log(100) to the logged outcome:
#       log(100 * wage) = log(wage) + log(100).
# So "relogging" the outcome in different units == adding the constant c=log(100).
# That is why the cents-vs-dollars comparison IS a translation-invariance test.
#
# ------------------------------------------------------------------------------
# METHOD B  — OUTCOME-WEIGHT (algebraic) check.   [does NOT rerun the estimator]
# ------------------------------------------------------------------------------
# Take the FIXED weight vector w from a single fit and apply it to both outcome
# codings. Report, to many decimals:
#   tau_dollars = sum(w * Y_dollars)
#   tau_cents   = sum(w * Y_cents)          (Y_cents = Y_dollars + k)
#   diff_actual = tau_cents - tau_dollars
#   diff_pred   = k * sum(w)                (the algebra's prediction)
#   sum_w       = sum(w)                    (= 0  iff translation invariant)
# Weights never change; only the outcome coding changes.
translation_check_weights <- function(name, w, Y_dollars, Y_cents, k,
                                       digits = 12) {
  tau_d <- sum(w * Y_dollars)
  tau_c <- sum(w * Y_cents)
  sumw  <- sum(w)
  data.frame(
    Estimator   = name,
    tau_dollars = round(tau_d, digits),
    tau_cents   = round(tau_c, digits),
    diff_actual = round(tau_c - tau_d, digits),
    diff_pred   = round(k * sumw,    digits),
    sum_w       = round(sumw,        digits),
    invariant   = isTRUE(all.equal(sumw, 0, tolerance = 1e-10)),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# METHOD A  — RERUN-THE-ESTIMATOR (behavioural) check.   [refits DML twice]
# ------------------------------------------------------------------------------
# Fit once on Y_dollars, fit AGAIN on the shifted/relogged outcome Y_cents,
# using the SAME seed both times so cross-fitting splits and forest randomness
# are identical (any difference is then a real invariance result, not RNG).
# Compare the realised shift tau(Y_cents) - tau(Y_dollars) against k * sum(w).
#
# fit_fun(Y) must return a list/obj from which `estimator_row` can be summarised;
# here we pass a closure that calls dml_with_smoother and returns summary() and
# the weights. Returns one row per estimator in `estimators`.
translation_check_rerun <- function(cell_name, fit_dollars, fit_cents,
                                     w_dollars, estimators, k, digits = 12) {
  res_d <- summary(fit_dollars)
  res_c <- summary(fit_cents)
  do.call(rbind, lapply(estimators, function(est) {
    tau_d <- get_estimate(res_d, est)
    tau_c <- get_estimate(res_c, est)
    w     <- w_dollars[[est]]
    data.frame(
      Cell        = cell_name,
      Estimator   = est,
      tau_dollars = round(tau_d, digits),
      tau_cents   = round(tau_c, digits),
      diff_actual = round(tau_c - tau_d, digits),  # realised shift from rerun
      diff_pred   = round(k * sum(w), digits),      # predicted shift = k*sum(w)
      sum_w       = round(sum(w), digits),
      invariant   = isTRUE(all.equal(tau_c, tau_d, tolerance = 1e-8)),
      stringsAsFactors = FALSE
    )
  }))
}


# ==============================================================================
# §K  UNIFIED TRANSLATION-INVARIANCE ROW (Method A + Method B in one row)
# ==============================================================================
# Produces a single table row carrying BOTH translation-invariance diagnostics
# for one estimator, so kappa / DML / learner estimators can be stacked in one
# consolidated table.
#
#   Method B (weight / algebraic):  diff_B = k * sum(w_dollars)
#       -> uses ONLY the dollars-fit weights; no refit. = 0  iff  sum(w)=0.
#   Method A (rerun / behavioural): diff_A = tau_cents - tau_dollars
#       -> the realised shift from estimating on both outcome codings.
#
# For weighted-sum estimators whose weights do NOT depend on Y (all kappa
# estimators), the two outcome codings reuse identical weights, so Method A and
# Method B coincide to machine precision. For adaptive learners (DML smoother,
# ranger/xgboost via DoubleML) the weights are re-learned, so diff_A may differ
# from diff_B; the gap (gap_AB) is the non-affine-smoother residual.
#
# Arguments:
#   name        : estimator label
#   w_dollars   : outcome-weight vector from the dollars fit
#   tau_dollars : point estimate on Y (dollars)
#   tau_cents   : point estimate on Y + k (cents)   [Method A rerun result]
#   k           : shift constant = log(100)
#   digits      : decimals to report (default 12)
ti_unified_row <- function(name, w_dollars, tau_dollars, tau_cents, k,
                           digits = 12) {
  sumw   <- sum(w_dollars)
  diff_B <- k * sumw                 # Method B prediction
  diff_A <- tau_cents - tau_dollars  # Method A realised shift
  data.frame(
    Estimator   = name,
    tau_dollars = round(tau_dollars, digits),
    tau_cents   = round(tau_cents,   digits),
    sum_w       = round(sumw,        digits),
    diff_B      = round(diff_B,      digits),  # Method B: k * sum(w)
    diff_A      = round(diff_A,      digits),  # Method A: rerun shift
    gap_AB      = round(diff_A - diff_B, digits),
    invariant_B = isTRUE(all.equal(sumw, 0, tolerance = 1e-10)),
    invariant_A = isTRUE(all.equal(tau_cents, tau_dollars, tolerance = 1e-8)),
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# END functions_all.R
# ==============================================================================
