################################################################################
##### My made of IV design 
###############################################################################
# Simulated IV data for debugging PLIV and IIVM

set.seed(123)

N <- 1000

X <- matrix(rnorm(N * 5), nrow = N, ncol = 5)
colnames(X) <- paste0("X", 1:5)

# Binary instrument
Z_prob <- plogis(0.5 * X[, 1] - 0.3 * X[, 2])
Z <- rbinom(N, 1, Z_prob)

# Binary treatment affected by Z
D_prob <- plogis(-0.2 + 1.0 * Z + 0.4 * X[, 1] - 0.2 * X[, 3])
D <- rbinom(N, 1, D_prob)

# Outcome affected by D and X
Y <- 2 * D + X[, 1] + 0.5 * X[, 2] + rnorm(N)

data_dml_iv <- double_ml_data_from_matrix(
  X = X,
  y = Y,
  d = D,
  z = Z
)
# this is my small IV design

# Now lets try to estimate it 
# First I activate both OLS (linear and logistic regression)
lrn_ols <- lrn("regr.lm")
lrn_logit <- lrn("classif.log_reg", predict_type = "prob")


# Now estimating the PLR IV using the instrument that I havent used before,
# I need to deine an additional nuisance paramter
pliv_parametric_dml <- DoubleMLPLIV$new(
  data_dml_iv,
  ml_l = lrn_ols,
  ml_m = lrn_ols,
  ml_r = lrn_ols,
  n_folds = 1,
  apply_cross_fitting = FALSE
)
# the functions wants an regression even though we have a binary variable.

# we fitt the modle
pliv_parametric_dml$fit(
  store_predictions = TRUE,
  store_models = TRUE
)

print(pliv_parametric_dml)

# lets see whetehr the outcome weights stuff works 
omega_pliv <- get_outcome_weights(pliv_parametric_dml, data_dml_iv)

all.equal(
  as.numeric(omega_pliv$omega %*% Y),
  as.numeric(pliv_parametric_dml$coef)
)




#####################################################################################
###### Lets check it for the Wald AIPW 
###############################################################################
# Wald-AIPW / IIVM with OLS outcome regression and logit nuisance learners

iivm_parametric_dml <- DoubleMLIIVM$new(
  data_dml_iv,
  ml_g = lrn("regr.lm"),
  ml_m = lrn("classif.log_reg", predict_type = "prob"),
  ml_r = lrn("classif.log_reg", predict_type = "prob"),
  n_folds = 1,
  apply_cross_fitting = FALSE
)

iivm_parametric_dml$fit(
  store_predictions = TRUE,
  store_models = TRUE
)

print(iivm_parametric_dml)


omega_iivm <- get_outcome_weights(iivm_parametric_dml, data_dml_iv)

all.equal(
  as.numeric(omega_iivm$omega %*% Y),
  as.numeric(iivm_parametric_dml$coef)
)



##### next test 
###############################################################################
# Unit test for get_smoother_weights_lm() inside the IIVM object: g1 model

# Take one of the fitted outcome models from IIVM
fold_models_g1 <- iivm_parametric_dml$models$ml_g1$d[[1]]
model_g1 <- fold_models_g1[[1]]$model

# Extract data
data_iv <- as.matrix(data_dml_iv$data)
x_cols <- data_dml_iv$x_cols
y_col <- data_dml_iv$y_col
Z <- data_iv[, data_dml_iv$z_cols]

# g1 is trained on observations with Z == 1
sub_z1 <- Z == 1

X_train <- data_iv[sub_z1, x_cols, drop = FALSE]
Y_train <- data_iv[sub_z1, y_col, drop = FALSE]
X_new   <- data_iv[, x_cols, drop = FALSE]

# Build smoother manually
S_g1_test <- get_smoother_weights_lm(
  model = model_g1,
  X = X_train,
  Y = Y_train,
  Xnew = X_new
)

# Use smoother to reconstruct predictions
Yhat_g1_from_S <- as.numeric(S_g1_test %*% as.numeric(Y_train))

# Compare to DoubleML prediction
all.equal(
  Yhat_g1_from_S,
  as.numeric(iivm_parametric_dml$predictions$ml_g1)
)


######## This could be good
###############################################################################
# Full test: custom outcome weights for Wald-AIPW / IIVM

omega_iivm_my <- my_get_outcome_weights_DoubleML(
  iivm_parametric_dml,
  data_dml_iv
)

all.equal(
  as.numeric(omega_iivm_my$omega %*% Y),
  as.numeric(iivm_parametric_dml$coef)
)






