rm(list=ls())
#check number twp
################################################################################
##### Simulated IV design
################################################################################
set.seed(123)
N <- 1000
X <- matrix(rnorm(N * 5), nrow = N, ncol = 5)
colnames(X) <- paste0("X", 1:5)

Z_prob <- plogis(0.5 * X[, 1] - 0.3 * X[, 2])
Z      <- rbinom(N, 1, Z_prob)

D_prob <- plogis(-0.2 + 1.0 * Z + 0.4 * X[, 1] - 0.2 * X[, 3])
D      <- rbinom(N, 1, D_prob)

Y <- 2 * D + X[, 1] + 0.5 * X[, 2] + rnorm(N)

data_dml_iv <- double_ml_data_from_matrix(X = X, y = Y, d = D, z = Z)

################################################################################
##### PLR-IV — OLS for all three nuisance params, n_folds = 5
################################################################################
pliv_cf <- DoubleMLPLIV$new(
  data_dml_iv,
  ml_l = lrn("regr.lm"),
  ml_m = lrn("regr.lm"),   # LPM for E[Z|X]
  ml_r = lrn("regr.lm"),   # LPM for E[D|X]
  n_folds = 5,
  apply_cross_fitting = TRUE
)
pliv_cf$fit(store_models = TRUE, store_predictions = TRUE)
print(pliv_cf)

# Test 1: published package
omega_pliv_pub <- get_outcome_weights(pliv_cf, data_dml_iv)
cat("PLR-IV published (n_folds=5):",
    all.equal(as.numeric(omega_pliv_pub$omega %*% Y),
              as.numeric(pliv_cf$coef)), "\n")

# Test 2: your custom function
omega_pliv_my <- my_get_outcome_weights_DoubleML(pliv_cf, data_dml_iv)
cat("PLR-IV my function (n_folds=5):",
    all.equal(as.numeric(omega_pliv_my$omega %*% Y),
              as.numeric(pliv_cf$coef)), "\n")

################################################################################
##### Wald-AIPW — OLS outcome + logit propensity, n_folds = 5
################################################################################
iivm_cf <- DoubleMLIIVM$new(
  data_dml_iv,
  ml_g = lrn("regr.lm"),
  ml_m = lrn("classif.log_reg", predict_type = "prob"),
  ml_r = lrn("classif.log_reg", predict_type = "prob"),
  n_folds = 5,
  apply_cross_fitting = TRUE
)
iivm_cf$fit(store_models = TRUE, store_predictions = TRUE)
print(iivm_cf)

# Test 1: published package
omega_iivm_pub <- get_outcome_weights(iivm_cf, data_dml_iv)
cat("Wald-AIPW published (n_folds=5):",
    all.equal(as.numeric(omega_iivm_pub$omega %*% Y),
              as.numeric(iivm_cf$coef)), "\n")

# Test 2: your custom function
omega_iivm_my <- my_get_outcome_weights_DoubleML(iivm_cf, data_dml_iv)
cat("Wald-AIPW my function (n_folds=5):",
    all.equal(as.numeric(omega_iivm_my$omega %*% Y),
              as.numeric(iivm_cf$coef)), "\n")

################################################################################
##### Unit test: get_smoother_weights_lm() on a fold-specific model
################################################################################
# Pull one fold-specific lm model from the n_folds=5 IIVM object
fold_models_g1 <- iivm_cf$models$ml_g1$d[[1]]
model_g1_raw   <- fold_models_g1[[1]]$model

# Unwrap mlr3 wrapping to get raw lm object
model_g1 <- if (inherits(model_g1_raw$model, "lm")) {
  model_g1_raw$model
} else {
  model_g1_raw
}

# Extract data
data_iv <- as.matrix(data_dml_iv$data)
x_cols  <- data_dml_iv$x_cols
y_col   <- data_dml_iv$y_col
Z_vec   <- data_iv[, data_dml_iv$z_cols]

# fold 1 test/train indices
test_ids  <- iivm_cf$smpls[[1]]$test_ids[[1]]
train_ids <- iivm_cf$smpls[[1]]$train_ids[[1]]

# g1 trained on Z==1 observations within training fold
sub_z1        <- Z_vec == 1
train_ids_sub <- train_ids[sub_z1[train_ids]]

X_train <- data_iv[train_ids_sub, x_cols, drop = FALSE]
Y_train <- data_iv[train_ids_sub, y_col,  drop = FALSE]
X_new   <- data_iv[test_ids,      x_cols, drop = FALSE]

# Build smoother for this fold
S_fold <- get_smoother_weights_lm(
  model = model_g1,
  X     = X_train,
  Y     = Y_train,
  Xnew  = X_new
)

# Reconstruct predictions from smoother
Yhat_from_S <- as.numeric(S_fold %*% as.numeric(Y_train))

# Compare to what DoubleML stored for these test observations
cat("Smoother unit test (fold 1, g1):",
    all.equal(Yhat_from_S,
              as.numeric(iivm_cf$predictions$ml_g1)[test_ids]), "\n")



methods("get_smoother_weights")
getAnywhere("get_smoother_weights.lm")
getAnywhere("get_smoother_weights.ranger")
getAnywhere("get_smoother_weights.xgb.Booster")
getAnywhere("get_smoother_weights.default")


