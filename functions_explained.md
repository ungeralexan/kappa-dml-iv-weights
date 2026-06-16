# Function Documentation — `functions_all.R`
For
### Thesis: *Comparing Kappa Weighting and Causal Machine Learning Estimators via Weight-Based Diagnostics*
 
> This document will show the central functions used in my Master Thesis and explain what it entails


 
---
 

## 1 `safe_logit()`
 
**Section:** §a | **Part:** Part 1 — Kappa Estimators
---
 
### The code
 
```r
safe_logit <- function(eta) {
  eta <- pmin(pmax(eta, -35), 35)
  1 / (1 + exp(-eta))
}
```

#### What it does 
`safe_logit()` takes a real number (or a vector of real numbers) and maps it to a probability between 0 and 1. It is the standard **logistic (sigmoid) function**. It adds one protection as it clips extreme values before computing, so the calculation never crashes due to floating-point overflow.
It is the **numerical primitive** shared by every propensity score routine in the file.
the input is a numeric scalar or a vector and the output are probability values between 0 and 1, numeric same size as the input eta.





## 2 `prep_design_for_mest()`

---
 
### The code
 
```r
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
```

#### What it does 
It takes the covariate matrix X and standardizes non-intercept column to have mean 0 and standard deviation 1.
The purpose is that when it comes to estimation we beed to compute the Jacobian and invert so we can get the standard errors, If teh Covariates have scales that are different that could make the Jacobian matrix ill conditioned. 
The inputs are a numeric matrix X or data frame. Must contain an intercept column (all 1s) plus covariate columns.
The output is a numeric matrix Same dimensions as `X`. All other columns transformed to mean 0, standard deviation 1.



## 3 `logit_mle()`

---
### The Code
```r
logit_mle <- function(Z, X) {
  df  <- data.frame(Z = Z, X)
  fit <- glm(Z ~ . - 1, data = df, family = binomial(link = "logit"))
  fitted.values(fit)
}

```
### What it does 
The input is Z which is a vector full of Z and 1 per obs, X which is the design covariate matrix. It already includes an intercept column (which is full of 1)
You receive a vector of probabilities, where each one stands for one observation.
They are between 0 and 1 as propensity scores. 
Important we surpress automatic the intercept. X should already have ones. The glm of R would already add another intercept thats why we surpress the automatic one. 
In Słoczyński, Uysal & Wooldridge (2025) mention that if X includes a constant the tau estimators all collapse to the same number, even though in final sample maybe not hold. 


## 4 `fit_logit_alpha()`

---
### The code
```r
fit_logit_alpha <- function(Z, X) {
  Z     <- as.numeric(Z)
  X     <- prep_design_for_mest(X)
  fit   <- glm.fit(x = X, y = Z, family = binomial(link = "logit"))
  alpha <- as.numeric(coef(fit));  alpha[is.na(alpha)] <- 0
  p     <- as.vector(pmin(pmax(safe_logit(X %*% alpha), 1e-8), 1 - 1e-8))
  list(alpha = alpha, p = p, X_used = X)
}
```

### What it does:
Gives coefficients, p and the standardizes X matrix. Will be needed for sandwhich SE formula.
It again takes the instrument vector, the covariate Matrix, inclusive intercept.
It outputs the estimated coeffs, the propensity score which is clipped and the standardised version of the X that was used to get the alphas. 
Important it returns X_used as the coefs (alpha) are on the standardized covariate sclae. If I want later compute X*alpha to reconstruct the linear predictor it is necessary to used the standardized not the origanl X. 
The clip  prevents p=0 or p=1 exactly, which would cause division by zero in the kappa formulas. 
We also prevent perfect collinearity and replace NA with zero which would mean that this varible has no effect. 


## 5 `fit_cbps_alpha()`

---
```r
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
```


### What it does:
This is a handcoded DBPS solver using Newtons method optimization loop.
We wanna find alha such that the wieghted covariate emans are balanced between Z=1 and Z=0.
Input are again instrument and and covariate matrix X a tolerance of 1e-9 and a maximum of iterations of 5000 same as in the SUW paper defined.
As output we have the CBPS coeffs on a standardized sclae, the CBPS propensity scores that are clipped. The X_used which is standardized. Moreover, it shows whether the algorithm converged and actually reached tolerance. Moreover the maximum absolut value of the moment conditions at the final iterate. Showing how close the exact balance is I got. 

The Newton Loop entails:
Initialising from the MLE solution as it is a already decnt guess. Will converge faster then starting from 0.
We first evaluate the moment condition and check how far are we from balance
If the norm is msaller than the best seen so far we save it as best_b
If the norm is below the tolerance we say convergence and stop
Then we compute the Jacobian (the matrix of derivatives of the moment conditions with respect to the coefficients) and solve J * step.
Then a backtracking line search is activated in which we try the full step then half, then quater up to 50 halvings until the momet norm actually decreases. 
I also implemented a ridge fallback as if the Jacobian is singular (collinearity), I add a small positive number along the diagonal (ridge) which makes the matrix invertible again



## 6 `cbps()`
---

```r
cbps <- function(Z, X, tol = 1e-9, max_iter = 5000, verbose = FALSE) {
  fit_cbps_alpha(Z, X, tol = tol, max_iter = max_iter)
}
```

### What it does : 
It is a wrapper, it just calls fit_cbps_alpha() and controls what gets handed back to the caller.
It exists for readability and inputs the instrument vector, the covariate mateix with intercpet and the parameters tol, max_iter.
It iotputs the same list from fit_cbps_alpha(): alpha, p, X_used, converged, max_moment



## 7 `get_cbps_p()`

---

```r
get_cbps_p <- function(Z, X) {
  out <- cbps(Z, X)
  if (is.list(out) && !is.null(out$p)) return(as.vector(out$p))
  as.vector(out)
}
```



### What it does 
As the next functions mainly need propensity scores and not alpha or the X matrix we just extract the out$p 
We again input the the instrument vector and teh covariate matrix with the intercepz and get as output the p vector whihc is a plain numeric vector of CBPS propensity scores. 


## 8 `kappa_outcome_weights()`

---

```r
kappa_weights <- function(Z, D, p) {
  list(
    kappa  = 1 - D * (1 - Z) / (1 - p) - (1 - D) * Z / p,
    kappa1 = D * (Z - p) / (p * (1 - p)),
    kappa0 = (1 - D) * ((1 - Z) - (1 - p)) / (p * (1 - p))
  )
}
```

### What it does:
Kappa ifentifies Abbadies original identification weight. Kappa 1 identifies the treated complier weight and kappa 0 identifies the untreated complier weight. 
It inputs the instruement (0/1), the treatment (0/1) and the proipensity scores
(maybe make sure they all are kind of binary).
It outputs a list of three vectors each length N :
- Kappa: WHich is abbadies original identification weight
- Kappa1: which is the treated complier weight
- Kappa0: which is the untreated complier weight


## 9 `tau_u()`

---

```r
tau_u <- function(Y, Z, D, p) {
  s1 <- sum(Z / p)
  s0 <- sum((1 - Z) / (1 - p))
  numerator   <- sum(Y * Z / p) / s1 - sum(Y * (1 - Z) / (1 - p)) / s0
  denominator <- sum(D * Z / p) / s1 - sum(D * (1 - Z) / (1 - p)) / s0
  numerator / denominator
}
```

### What it does:
This is the normalized estimator by UYsal which is translation invariance at the same time.
It Computes separately normalised IPW means for Z=1 and Z=0, then takes.
It takes as input the the outcome vector Y the unstrument the treatement and the propesnity scores. 
It gives as uotput a single number which is the estimator. 


## 10 `tau_a10()`

---

```r
tau_a10 <- function(Y, Z, D, p) {
  kw <- kappa_weights(Z, D, p)
  sum(kw$kappa1 * Y) / sum(kw$kappa1) - sum(kw$kappa0 * Y) / sum(kw$kappa0)
}
```

### What it does
It depicts : BADIE-CATTANEO NORMALIZED ESTIMATOR  [translation invariant]. It sepperately normalizes kappa 1 and kappa 0 weighted outcome weights. 
It takes the same inputs as the function above and outputs the LE using the Badie Cattaneo normalization



## 11 `tau_unnorm()`

---

```r
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
```


### What it does :
We use the property that all thre unnormlised estimators share the same numerator only the denominator differs. So rather than writing three separate functions, your code implements all three in one, with a switch statement to pick the denominator.
The input are the outcome, the instrument the treatment and the propensity score and a which function whihc determines which of the three functions I should use. Like which of three denominators and then gives as output the estimator we select.


## 12 `kappa_outcome_weights()`

---

```r
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
```

### What it does:
Every single estimator can be written as a sum of the outcomes times the weighst per observation. 
We take the five estimator formulas you already know and algebraically collapse each one into its weight form. 
The inputs are the Instrument, the treatment and he propensity score and as outcomes we get different weight vectors of length N.
Which is the number of observations
