library(DoubleML)
library(mlr3)
library(mlr3learners)
library(hdm)
library(OutcomeWeights)

data(pension)

Y <- pension$net_tfa
D <- pension$p401
X <- model.matrix(~ 0 + age + db + educ + fsize + hown + inc + male + marr + pira + twoearn,
                  data = pension)

data_dml <- double_ml_data_from_matrix(X = X, y = Y, d = D)

lrn_ols <- lrn("regr.lm")
lrn_logit <- lrn("classif.log_reg", predict_type = "prob")

plr_logit_dml <- DoubleMLPLR$new(
  data_dml,
  ml_l = lrn_ols,
  ml_m = lrn_logit,
  n_folds = 1,
  apply_cross_fitting = FALSE
)

plr_logit_dml$fit(store_predictions = TRUE, store_models = TRUE)

print(plr_logit_dml)

omega_plr <- get_outcome_weights(plr_logit_dml, data_dml)

all.equal(
  as.numeric(omega_plr$omega %*% Y),
  as.numeric(plr_logit_dml$coef)
)


###############################################################################
# Specification for the Wald AIPW
aipw_parametric_dml <- DoubleMLIRM$new(
  data_dml,
  ml_g = lrn_ols,
  ml_m = lrn_logit,
  n_folds = 1,
  apply_cross_fitting = FALSE
)

aipw_parametric_dml$fit(store_predictions = TRUE, store_models = TRUE)

omega_aipw <- get_outcome_weights(aipw_parametric_dml, data_dml)

all.equal(
  as.numeric(omega_aipw$omega %*% Y),
  as.numeric(aipw_parametric_dml$coef)
)


###############################################################################
##### checking whetehr get_double_mlsmoother works with reg.lm
fold_models <- plr_logit_dml$models$ml_l$d[[1]]
model <- fold_models[[1]]$model

class(model)
class(model$model)
str(model, max.level = 3)

# the covariate order to be keep in mind
lm_x_cols <- attr(model$terms, "term.labels")


###############################################################################
# Check whether lm fitted values equal DoubleML ml_l predictions

fold_models <- plr_logit_dml$models$ml_l$d[[1]]
model <- fold_models[[1]]$model

all.equal(
  as.numeric(model$fitted.values),
  as.numeric(plr_logit_dml$predictions$ml_l)
)


###############################################################################
# Check covariate ordering inside lm

attr(model$terms, "term.labels") # this the order used in the lm formula 
data_dml$x_cols # this is the order stored in the Double ML data object
# The lm model internally uses a different column order than data_dml$x_cols.

# in order to reconstruct the lm predictions I guess the safer source is : 
attr(model$terms, "term.labels")
# as this tells me what the fitted model actually used.


###############################################################################
# Manually construct OLS smoother and check whether S %*% Y equals fitted values

lm_x_cols <- attr(model$terms, "term.labels") # the exact variable order extracted used by the fitted lm model

X_lm <- as.matrix(data_dml$data[, ..lm_x_cols]) # builds covariate matrix in the same order as the lm model
Y_lm <- as.numeric(data_dml$data[[data_dml$y_col]]) # extracts the outcome verctor
X_aug <- cbind("(Intercept)" = 1, X_lm) # now adding the intercept

S_lm <- X_aug %*% solve(crossprod(X_aug)) %*% t(X_aug) # compute : S = X(X'X)^(-1)X'

Yhat_from_S <- as.numeric(S_lm %*% Y_lm) # this computes : Y_hat = S Y

all.equal(
  Yhat_from_S,
  as.numeric(model$fitted.values)
)




###############################################################################
# Same smoother using data_dml$x_cols order

dml_df <- as.data.frame(data_dml$data)

X_wrong <- as.matrix(dml_df[, data_dml$x_cols, drop = FALSE])
X_wrong_aug <- cbind("(Intercept)" = 1, X_wrong)

S_wrong <- X_wrong_aug %*% solve(crossprod(X_wrong_aug)) %*% t(X_wrong_aug)

Yhat_wrong <- as.numeric(S_wrong %*% Y_lm)

all.equal(
  Yhat_wrong,
  as.numeric(model$fitted.values)
)
# still works (somehow??)



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






###### Lets try my function

###############################################################################
# Unit test for get_smoother_weights_lm() inside the IIVM object

# Take one of the fitted outcome models from IIVM
fold_models_g0 <- iivm_parametric_dml$models$ml_g0$d[[1]]
model_g0 <- fold_models_g0[[1]]$model

# Extract data
data_iv <- as.matrix(data_dml_iv$data)
x_cols <- data_dml_iv$x_cols
y_col <- data_dml_iv$y_col
Z <- data_iv[, data_dml_iv$z_cols]

# g0 is trained on observations with Z == 0
sub_z0 <- Z == 0

X_train <- data_iv[sub_z0, x_cols, drop = FALSE]
Y_train <- data_iv[sub_z0, y_col, drop = FALSE]
X_new   <- data_iv[, x_cols, drop = FALSE]

# Build smoother manually
S_g0_test <- get_smoother_weights_lm(
  model = model_g0,
  X = X_train,
  Y = Y_train,
  Xnew = X_new
)

# Use smoother to reconstruct predictions
Yhat_g0_from_S <- as.numeric(S_g0_test %*% as.numeric(Y_train))

# Compare to DoubleML prediction
all.equal(
  Yhat_g0_from_S,
  as.numeric(iivm_parametric_dml$predictions$ml_g0)
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










###### try it for my example 

data(pension) # Find variable description if you type ?pension in console

# Treatment
D = pension$p401
# Instrument
Z = pension$e401
# Outcome
Y = pension$net_tfa
# Controls
X = model.matrix(~ 0 + age + db + educ + fsize + hown + inc + male + marr + pira + twoearn, data = pension)
var_nm = c("Age","Benefit pension","Education","Family size","Home owner","Income","Male","Married","IRA","Two earners")
colnames(X) = var_nm
data_dml <- double_ml_data_from_matrix(X = X, y = Y, d = D, z = Z)


iivm_parametric_dml <- DoubleMLIIVM$new(
  data_dml,
  ml_g = lrn("regr.lm"),
  ml_m = lrn("classif.log_reg", predict_type = "prob"),
  ml_r = lrn("classif.log_reg", predict_type = "prob"),
  n_folds = 1,
  apply_cross_fitting = FALSE
)


iivm_parametric_dml$fit(store_predictions = TRUE, store_models = TRUE)

print(iivm_parametric_dml)

omega_plr <- get_outcome_weights(iivm_parametric_dml, data_dml)

all.equal(
  as.numeric(omega_plr$omega %*% Y),
  as.numeric(iivm_parametric_dml$coef)
)



