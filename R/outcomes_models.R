#' @title Estimate Outcome Models via K-Fold Cross-Fitting
#'
#' @description
#' This function implements K-fold cross-fitting to estimate conditional
#' expectations for both egos and alters. The splitting is done at the
#' ego-network level to ensure that egos and their corresponding alters are
#' always in the same fold (either training or hold-out).
#'
#' This function estimates:
#' 1. Alter Model: E[Y | F, X, i in R_a] -> mu_a_1, mu_a_0
#' 2. Ego Model:   E[Y | Z, X, i in R_e] -> mu_e_1, mu_e_0
#'
#' @param Y_e A numeric vector of outcomes for the egos (length n_e).
#' @param Y_a A numeric vector of outcomes for the alters (length n_a).
#' @param X_e A numeric matrix of covariates for the egos (n_e rows).
#'            Can be `NULL` if no covariates are used.
#' @param X_a A numeric matrix of covariates for the alters (n_a rows).
#'            Must have the same number of columns as `X_e`. Can be `NULL`.
#' @param Z_e A binary numeric vector (0 or 1) of treatment assignments for egos.
#' @param F_a A binary numeric vector (0 or 1) of `observed` exposures for alters.
#' @param ego_id_a A numeric vector mapping each alter to their ego's *index*
#'        (an integer from 1 to n_e).
#' @param reg_model_egos A function for the egos' outcome model (e.g., `glm`).
#' @param reg_model_alters A function for the alters' outcome model (e.g., `glm`).
#' @param formula_egos A formula for the egos' model. Vars: `Y`, `Z`, `X1`, `X2`...
#' @param formula_alters A formula for the alters' model. Vars: `Y`, `F`, `X1`, `X2`...
#' @param n_folds The number of folds for cross-fitting (e.g., 2).
#' @param ... Additional arguments passed to the regression model functions
#'        (e.g., `family = binomial`).
#'
#' @return A list containing four numeric vectors of cross-fitted predictions:
#' \describe{
#'   \item{mu_a_1}{Predictions for alters with F=1 (length n_a).}
#'   \item{mu_a_0}{Predictions for alters with F=0 (length n_a).}
#'   \item{mu_e_1}{Predictions for egos with Z=1 (length n_e).}
#'   \item{mu_e_0}{Predictions for egos with Z=0 (length n_e).}
#' }
#' @keywords internal
#'
estimate_outcome_models_cf <- function(Y_e,
                                       Y_a,
                                       X_e = NULL,
                                       X_a = NULL,
                                       Z_e,
                                       F_a,
                                       ego_id_a,
                                       reg_model_egos,
                                       reg_model_alters,
                                       formula_egos,
                                       formula_alters,
                                       n_folds = 2,
                                       ...) {

  # --- Input Validation ---
  n_e <- length(Y_e)
  n_a <- length(Y_a)

  if (is.null(reg_model_egos) || !is.function(reg_model_egos)
      || is.null(reg_model_alters) || !is.function(reg_model_alters)
      || is.null(formula_egos) || !inherits(formula_egos, "formula")
      || is.null(formula_alters) || !inherits(formula_alters, "formula")) {
    warning("Invalid regression model functions or formulas. Setting predicted values to zero")
    return(list(
      mu_a_1 = numeric(n_a),
      mu_a_0 = numeric(n_a),
      mu_e_1 = numeric(n_e),
      mu_e_0 = numeric(n_e)
    ))
  }

  if (n_e == 0 || n_a == 0) {
    stop("Y_e and Y_a must not be empty.")
  }
  if (length(ego_id_a) != n_a) {
    stop("ego_id_a length must match Y_a length (n_a).")
  }
  if (!is.null(X_e) && (nrow(X_e) != n_e)) {
    stop("X_e must have n_e rows.")
  }
  if (!is.null(X_a) && (nrow(X_a) != n_a)) {
    stop("X_a must have n_a rows.")
  }
  if (!is.null(X_e) && !is.null(X_a) && (ncol(X_e) != ncol(X_a))) {
    stop("X_e and X_a must have the same number of columns.")
  }

  # --- Helper function to prepare data.frames for modeling ---
  # This handles NULL X and assigns generic col names for formulas
  prepare_model_data <- function(Y, Tr, X, tr_name) {
    dat <- data.frame(Y = Y)
    dat[[tr_name]] <- Tr

    if (!is.null(X)) {
      X_df <- as.data.frame(X)
      p <- ncol(X_df)
      if (p > 0) {
        colnames(X_df) <- paste0("X", 1:p)
        dat <- cbind(dat, X_df)
      }
    }
    return(dat)
  }

  # Helper to prepare "newdata" for prediction
  prepare_newdata <- function(X, tr_val, tr_name) {
    if (is.null(X)) {
      dat <- data.frame(tr_val)
      colnames(dat) <- tr_name
    } else {
      X_df <- as.data.frame(X)
      p <- ncol(X_df)
      if (p > 0) {
        colnames(X_df) <- paste0("X", 1:p)
      }
      dat <- cbind(tr_val, X_df)
      colnames(dat)[1] <- tr_name
    }
    return(dat)
  }


  # --- Initialize Output Vectors ---
  # These will be filled in the original order
  mu_a_1 <- numeric(n_a)
  mu_a_0 <- numeric(n_a)
  mu_e_1 <- numeric(n_e)
  mu_e_0 <- numeric(n_e)

  # --- Create Folds (at ego level) ---
  # Ensure folds are roughly equal and randomized
  fold_ids_e <- sample(cut(seq(1, n_e), breaks = n_folds, labels = FALSE))

  # Assign alters to the same fold as their ego
  fold_ids_a <- fold_ids_e[ego_id_a]


  # --- K-Fold Cross-Fitting Loop ---
  for (k in 1:n_folds) {
    # 1. Define training and hold-out indices
    idx_e_holdout <- which(fold_ids_e == k)
    idx_e_train <- which(fold_ids_e != k)
    idx_a_holdout <- which(fold_ids_a == k)
    idx_a_train <- which(fold_ids_a != k)

    # Check for empty folds which can happen with small n
    if (length(idx_e_train) == 0 ||
        length(idx_a_train) == 0 ||
        length(idx_e_holdout) == 0) {
      warning(paste(
        "Fold", k, "is empty or results in empty training/holdout set. Skipping.",
        "Consider using fewer folds or stratified sampling if n is small."
      ))
      next
    }

    # 2. Prepare Training Data
    dat_e_train <-
      prepare_model_data(Y_e[idx_e_train], Z_e[idx_e_train],
                         if (is.null(X_e))
                           NULL
                         else
                           X_e[idx_e_train, , drop = FALSE], "Z")

    dat_a_train <-
      prepare_model_data(Y_a[idx_a_train], F_a[idx_a_train],
                         if (is.null(X_a))
                           NULL
                         else
                           X_a[idx_a_train, , drop = FALSE], "F")

    # 3. Fit Models on Training Data
    fit_e <-
      do.call(reg_model_egos,
              list(formula = formula_egos, data = dat_e_train, ...))
    fit_a <-
      do.call(reg_model_alters,
              list(formula = formula_alters, data = dat_a_train, ...))


    # 4. Prepare Hold-out Data for Prediction

    # --- For Ego Model (mu_e_1, mu_e_0) ---
    X_e_holdout <-
      if (is.null(X_e))
        NULL
    else
      X_e[idx_e_holdout, , drop = FALSE]
    newdata_e_z1 <- prepare_newdata(X_e_holdout, 1, "Z")
    newdata_e_z0 <- prepare_newdata(X_e_holdout, 0, "Z")

    # --- For Alter Model (mu_a_1, mu_a_0) ---
    X_a_holdout <-
      if (is.null(X_a))
        NULL
    else
      X_a[idx_a_holdout, , drop = FALSE]
    newdata_a_f1 <- prepare_newdata(X_a_holdout, 1, "F")
    newdata_a_f0 <- prepare_newdata(X_a_holdout, 0, "F")


    # 5. Predict on Hold-out and Store in Final Vectors

    # Store predictions in the correct (original) index locations
    if (length(idx_e_holdout) > 0) {
      mu_e_1[idx_e_holdout] <-
        predict(fit_e, newdata = newdata_e_z1, type = "response")
      mu_e_0[idx_e_holdout] <-
        predict(fit_e, newdata = newdata_e_z0, type = "response")
    }

    if (length(idx_a_holdout) > 0) {
      mu_a_1[idx_a_holdout] <-
        predict(fit_a, newdata = newdata_a_f1, type = "response")
      mu_a_0[idx_a_holdout] <-
        predict(fit_a, newdata = newdata_a_f0, type = "response")
    }
  }

  # --- 6. Return Results ---
  results <- list(
    mu_a_1 = mu_a_1,
    mu_a_0 = mu_a_0,
    mu_e_1 = mu_e_1,
    mu_e_0 = mu_e_0
  )

  # Final check for any NA/NaN (e.g., from failed folds)
  if (any(sapply(results, function(v)
    any(!is.finite(v))))) {
    warning(
      "Cross-fitting produced non-finite (NA/NaN) predictions. ",
      "This may be due to small n, empty folds, or model instability."
    )
  }

  return(results)
}
