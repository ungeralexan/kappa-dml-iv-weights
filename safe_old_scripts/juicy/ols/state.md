# Current state of the implementation

I tried to reproduce the Introduction to Causal ML estimators in R Notebook, using the instrumental varible versions, DoubleMLIIVM$new for the Wald AIPW and DoubleMLPLIV$new for the PLR- IV.
When trying to etsimate the nuisance learners especially E[D|X] and E[Z|X] in the PLR-Iv uisng logistic regression as the tretament and instrument are bianry I got an error. 
The reason is that DoubleMLPLIV internally enforces that all three nuisance learners (ml_l, ml_m, ml_r) must be regression tasks. 
The package does not support classifier for PLR IV even if Z is binary.
On the other hand DoubleMLIIVM, I can use the logit for the binary nuisance learners ml_m and ml_r because IIVM is designed for bianry Z and D and the package permits it.

## What could be my choices?
For the PLR-IV use regr.lm for everything. Then this would include treating E[Z|X] as a linear probability model and treat E[D|X] as well as a linear probability model.
Or a 2 stage workaround? must be fitting logit outside DoubleML for ml_m and ml_r, get the residuals, then feed the residualized Z and D into a simpler DoubleML structure. -> Would be quite complex.


## The smoother function calculation for R
```r
get_smoother_weights_lm <- function(model, X, Y, Xnew) {
  X = as.data.frame(X)
  Xnew = as.data.frame(Xnew)
  
  vars = attr(model$terms, "term.labels")
  X = X[, vars, drop = FALSE]
  Xnew = Xnew[, vars, drop = FALSE]
  
  has_intercept = attr(model$terms, "intercept") == 1
  X_mat = as.matrix(X)
  Xnew_mat = as.matrix(Xnew)
  
  if (has_intercept) {
    X_mat = cbind("(Intercept)" = 1, X_mat)
    Xnew_mat = cbind("(Intercept)" = 1, Xnew_mat)
  }
  
  S = Xnew_mat %*% MASS::ginv(crossprod(X_mat)) %*% t(X_mat) # this is the soother matrix used in the fomrmula
  return(S)
}
```

#### A deeper insight into the intercept
When lrn("regr.lm") fits inside DoubleML, the data is passed as a plain matrix with column names like X1, X2,.. There is no intercept column. The lm() call inside the mlr3 calls the inetrcept automatically via the formula y ~  which will means the condition will be true and the code will preprend a column of ones. 

In my framework X and Xnew coming into get_smoother_weights_lm are plain covariates matrices without intercoet. The intecept is then added isnide the function

Important Y is not used at all we dont need it, but I keep it to match the API of get_smoother_weights()


### What about my personal get smoother function
I build an NxN smooter matrix S by looping over cross fitting folds. For each fold I compzte a fold-specific smoother block and places it in the right position of S. This function is only called fro the outcome model only

```r
my_get_DoubleML_smoother <- function(fold_models,
                                  test_ids_list,
                                  train_ids_list,
                                  data,
                                  x_cols,
                                  y_col,
                                  subset = NULL) {
  N <- nrow(data)
  S <- matrix(0, N, N)
  
  if (is.null(subset)) subset <- rep(TRUE, N)
  
  for (i in seq_along(fold_models)) {
    model <- fold_models[[i]]$model
    
    # ranger's double-wrapping quirk: 
    if (inherits(model$model, "ranger")) model <- model$model
    
    # xgboost's feature ordering quirk: 
    if (inherits(model, "xgb.Booster")) x_cols <- xgboost::getinfo(model, "feature_name")
    
    
    test_ids  <- test_ids_list[[i]]
    train_ids <- train_ids_list[[i]]
    train_ids_sub <- train_ids[subset[train_ids]]
    
    test_data   <- data[test_ids,      x_cols, drop = FALSE]
    train_data  <- data[train_ids_sub, x_cols, drop = FALSE]
    train_label <- data[train_ids_sub, y_col,  drop = FALSE]
    
    if (inherits(model, "lm")) {
      S_fold <- get_smoother_weights_lm(
        model = model,
        X = train_data,
        Y = train_label,
        Xnew = test_data
      )
    } else {
      S_fold <- OutcomeWeights:::get_smoother_weights(
        model,
        X = train_data,
        Y = train_label,
        Xnew = test_data
      )
    }
    
    S[test_ids, train_ids_sub] <- S_fold
  }
  
  return(S)
}
```
Why do ranger and XGBOOST checks come before data extraction.
Ranger modfies model itself  it unwraps a double-wrapped object. If I would extrat data first and tehn unwrap the model passed downstream would be wrong.
For XGBoost it modifies x_cols which must happen before data[test_ids, x_cols, drop=FALSE] otherwhis I would be selecting columns with the wrong names and crash or get wrong data.

LM does not need to modify or unwrappe anything.
Important y_col is used to extract the training outcomes for the fold. This gets passed as Y to get_smoother_weights_lm.
Interesting for OLS Y not used in computing the matrix as we saw but by ranger and XGboost is used inside OutcomeWeights:::get_smoother_weights, because those tree-based methods need the outcomes to compute their smoother weights.

Subset is also important. It is a logical vector of length N (e.g. Z == 0 or Z == 1). The line train_ids[subset[train_ids]] filters the training indices to keep only those where the subset condition is TRUE — i.e. only the Z=0 (or Z=1) training observations.
train_data only conatins Z = 0 or Z = 1. The smoother SfoldS_{fold}
Sfold​ maps from those subset training observations to all test observations
The reason for that is that for the Wald AIPW, outcome models are only trained on one side of the instrument :
- ml_g0 = E^[Y∣Z=0,X]\hat{E}[Y|Z=0, X]
E^[Y∣Z=0,X] → trained only on observations where Z=0Z=0
Z=0
- ml_g1 = E^[Y∣Z=1,X]\hat{E}[Y|Z=1, X]
E^[Y∣Z=1,X] → trained only on observations where Z=1Z=1
Z=1


## The get outcome weights function using the S matrix
```r
my_get_outcome_weights_DoubleML = function(object, dml_data, ...) {
  ## Preps
  type  <- class(object)[1] # which estimator "DoubleMLPLIV", "DoubleMLIIVM"
  data  <- as.matrix(dml_data$data)
  preds <- lapply(object$predictions, as.numeric) # here the nuisance pred stored by double ml 
  
  y_col  <- dml_data$y_col
  x_cols <- dml_data$x_cols
  d_col  <- dml_data$d_cols
  z_col  <- dml_data$z_cols
  
  Y <- data[, y_col]
  X <- data[, x_cols, drop = FALSE]
  D <- if (!is.null(d_col)) data[, d_col]
  Z <- if (!is.null(z_col)) data[, z_col]
  N <- nrow(data)
  
  test_ids_list  <- object$smpls[[1]]$test_ids # list of K vectors 
  train_ids_list <- object$smpls[[1]]$train_ids
  
  
  ## Get outcome weights
  ### PLR ###
  if (type == "DoubleMLPLR") {
    fold_models <- object$models$ml_l$d[[1]]
    
    S <- my_get_DoubleML_smoother(
      fold_models, test_ids_list, train_ids_list, data, x_cols, y_col, NULL
    )
    
    Z.tilde = D.tilde = D - preds$ml_m
    T_mat = diag(N) - S
    omega = OutcomeWeights:::pive_weight_maker(Z.tilde, D.tilde, T_mat)
    
    
    ### PLR-IV ###
  } else if (type == "DoubleMLPLIV") {
    fold_models <- object$models$ml_l$d[[1]] # outcome model E[Y|X]
    
    S <- my_get_DoubleML_smoother(
      fold_models, test_ids_list, train_ids_list, data, x_cols, y_col, NULL # getting the smoother for ml_l only 
    )
    
    Z.tilde = Z - preds$ml_m  # instrument residual: Z̃ = Z - Ê[Z|X]
    D.tilde = D - preds$ml_r # treatment residual:  D̃ = D - Ê[D|X
    T_mat = diag(N) - S
    omega = OutcomeWeights:::pive_weight_maker(Z.tilde, D.tilde, T_mat)
    
    
    ### AIPW-ATE ###
  } else if (type == "DoubleMLIRM") {
    fold_models_g0 <- object$models$ml_g0$d[[1]]
    fold_models_g1 <- object$models$ml_g1$d[[1]]
    
    sub_d0 <- D == 0
    sub_d1 <- D == 1
    
    S.d0 <- my_get_DoubleML_smoother(
      fold_models_g0, test_ids_list, train_ids_list, data, x_cols, y_col, sub_d0
    )
    S.d1 <- my_get_DoubleML_smoother(
      fold_models_g1, test_ids_list, train_ids_list, data, x_cols, y_col, sub_d1
    )
    
    Z.tilde = D.tilde = rep(1, N)
    lambda1 = D / preds$ml_m
    lambda0 = (1 - D) / (1 - preds$ml_m)
    T_mat = S.d1 - S.d0 + 
      lambda1 * (diag(N) - S.d1) - 
      lambda0 * (diag(N) - S.d0)
    omega = OutcomeWeights:::pive_weight_maker(Z.tilde, D.tilde, T_mat)
    
    
    ### Wald-AIPW ###
  } else if (type == "DoubleMLIIVM") { # two smoothers, one per instrument value
    fold_models_g0 <- object$models$ml_g0$d[[1]] # E[Y|Z=0,X] smoother
    fold_models_g1 <- object$models$ml_g1$d[[1]] # E[Y|Z=1,X] smoother
    
    sub_z0 <- Z == 0
    sub_z1 <- Z == 1
    
    S.z0 <- my_get_DoubleML_smoother(
      fold_models_g0, test_ids_list, train_ids_list, data, x_cols, y_col, sub_z0
    )
    S.z1 <- my_get_DoubleML_smoother(
      fold_models_g1, test_ids_list, train_ids_list, data, x_cols, y_col, sub_z1
    )
    
    lambdaz1 = Z / preds$ml_m
    lambdaz0 = (1 - Z) / (1 - preds$ml_m)
    Z.tilde = rep(1, N)
    D.tilde = preds$ml_r1 - preds$ml_r0 +  # Ê[D|Z=1,X] - Ê[D|Z=0,X]
      lambdaz1 * (D - preds$ml_r1) -   # IPW correction Z=1 side
      lambdaz0 * (D - preds$ml_r0) # IPW correction Z=0 side
    T_mat = S.z1 - S.z0 + 
      lambdaz1 * (diag(N) - S.z1) -  # Z=1 arm residualizer
      lambdaz0 * (diag(N) - S.z0)    # Z = 0 arm residualizer
    omega = OutcomeWeights:::pive_weight_maker(Z.tilde, D.tilde, T_mat)
    
  } else {
    stop(
      "Outcome weights are extracted for DoubleMLPLR, DoubleMLPLIV, DoubleMLIIVM, and DoubleMLIRM models"
    )
  }
  
  output = list("omega" = omega, "treat" = D)
  class(output) = c("get_outcome_weights")
  return(output)
}
```



## The distinction.
World one - > dml_with_smoother()
This is the simplest. Knaus wrote his own DML from scratch. You give it raw vectors:
```r
dml_with_smoother(Y, D, X, Z, n_cf_folds = 5)
```
Internally it does everything itself — fits ranger, stores the smoother matrix S as it goes, computes the weights. You never see any of the internals. It spits out one object and you call get_outcome_weights() on it. Done.
Limitation it only uses rnager, no OLS, no logit and hence is hardcoded



World 2 — > His DML_smoothers__2_.rmd notebook — the DoubleML bridge
 Here he uses the DoubleML package (separate package, written by other people). This package is much more general — you can plug in any learner.
 The mysterious functions you mentioned — make_plr_CCDDHNR2018, make_pliv_CHS2015, make_iivm_data — these are just fake data generators for testing. CCDDHNR are the authors' initials (Chernozhukov, Chetverikov, Demirer, Duflo, Hansen, Newey, Robins). CHS = Chernozhukov, Hansen, Spindler. He uses fake data just to demonstrate his package works. Nothing to do with Vietnam or 401k.


Step 1: create fake test data
        make_pliv_CHS2015(1000, dim_x=10)  ← just N=1000 fake observations

Step 2: define learners
        lrn("regr.lm")   ← OLS

Step 3: create DoubleML object
        DoubleMLPLIV$new(pliv_data, ml_l, ml_m, ml_r, n_folds=3)

Step 4: fit it
        pliv_obj$fit(store_models = TRUE, store_predictions = TRUE)
        ↑ this is CRITICAL — store_models saves each fold's lm object

Step 5: extract weights using his bridge function
        get_outcome_weights(object = pliv_obj, dml_data = pliv_data)
        ↑ this calls get_DoubleML_smoother() internally
        ↑ which loops over stored fold models and rebuilds S post-hoc

Step 6: check
        omega %*% Y == coef  →  TRUE?


The thing is that get_DoubleML_smoother() in the published package doesn't handle lm properly with real cross-fitting. It works in his notebook because he might be using the dev version, or because certain conditions make it accidentally work.


World 3 — Your implementation — extending the bridge
Your world is identical to World 2 in structure, but you replace his get_DoubleML_smoother() with your my_get_DoubleML_smoother() which explicitly handles lm. And crucially you use it on real data (Vietnam), not fake test data. Your flow:

Step 1: create DoubleMLData from real Vietnam data
        DoubleMLData$new(df, y_col="y", d_cols="d", z_cols="z", ...)

Step 2: define learners — THE NEW PART
        lrn("regr.lm")                              ← OLS for outcome
        lrn("classif.log_reg", predict_type="prob") ← logit for propensity

Step 3: create DoubleML object — same as World 2
        DoubleMLPLIV$new(dat_cub, ml_l, ml_m, ml_r, n_folds=5)
        DoubleMLIIVM$new(dat_cub, ml_g, ml_m, ml_r, n_folds=5)

Step 4: fit — same as World 2
        obj$fit(store_models = TRUE, store_predictions = TRUE)

Step 5: extract weights — YOUR function instead of his
        my_get_outcome_weights_DoubleML(obj, dat_cub)
        ↑ calls my_get_DoubleML_smoother()
        ↑ which hits your new if(inherits(model, "lm")) branch
        ↑ which calls get_smoother_weights_lm()
        ↑ which computes Xnew(X'X)^{-1}X'  correctly

Step 6: check — same goal
        omega %*% Y == coef  →  TRUE  ← this is the proof

