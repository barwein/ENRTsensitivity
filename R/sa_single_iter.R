# source("src/bias_adjustment.r")
# source("src/outcomes_models.R")

#' @title Run a Single Iteration of the Sensitivity Analysis
#'
#' @description
#'  This function performs one iteration of the sensitivity analysis for a given
#' dataset (or bootstrap sample). It first estimates the outcome models using
#' cross-fitting and then uses those predictions to calculate the augmented,
#' bias-corrected estimates for IE and DE across a grid of sensitivity parameters.
#'
#' @param Y_e A numeric vector of outcomes for the egos.
#' @param Y_a A numeric vector of outcomes for the alters.
#' @param X_e A numeric matrix of covariates for the egos.
#' @param X_a A numeric matrix of covariates for the alters.
#' @param Z_e A binary numeric vector (0 or 1) of treatment assignments for egos.
#' @param F_a A binary numeric vector (0 or 1) of `observed` exposures for alters.
#' @param ego_id_a A numeric vector mapping each alter to their ego's index
#'        (1 to n_e). Required for IE variance.
#' @param reg_model_egos A function for the egos' outcome regression model.
#' @param reg_model_alters A function for the alters' outcome regression model.
#' @param formula_egos A formula for the egos' regression model.
#' @param formula_alters A formula for the alters' regression model.
#' @param pi_list_ego_ego A list where each element is a vector of exposure
#'        probabilities (`pi_i^e`) for egos.
#' @param pi_list_alter_ego A list where each element is a vector of exposure
#'        probabilities (`pi_i^a`) for alters.
#' @param kappa_vec A numeric vector of values for the kappa sensitivity parameter.
#' @param pz The known probability of an ego being assigned to treatment, Pr(Z=1).
#' @param n_folds The number of folds for cross-fitting outcome models.
#' @param estimate_var Logical. If TRUE, compute and return empirical variance
#'        estimates. If FALSE, variance columns will be NA.
#' @param ... Additional arguments passed to the regression model functions.
#'
#' @return A list containing two data.tables:
#' \describe{
#'   \item{IE_corrected}{A data.table with columns: `pi_param`, `ie_rd`, `ie_rd_var`}
#'   \item{DE_corrected}{A data.table with columns: `pi_param`, `kappa`, `de_rd`, `de_rd_var`}
#' }
#'
#' @keywords internal
#'
SA_one_iter <- function(Y_e,
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
                        pi_list_ego_ego,
                        pi_list_alter_ego,
                        kappa_vec,
                        pz,
                        n_folds = 2,
                        estimate_var = FALSE,
                        ...){

  # Check inputs
  if (!is.vector(Y_e) || !is.numeric(Y_e)) stop("Y_e must be a numeric vector.")
  if (!is.vector(Y_a) || !is.numeric(Y_a)) stop("Y_a must be a numeric vector.")
  if (!is.null(X_e) && (!is.matrix(X_e) || !is.numeric(X_e))) stop("X_e must be a numeric matrix.")
  if (!is.null(X_e) && (!is.matrix(X_a) || !is.numeric(X_a))) stop("X_a must be a numeric matrix.")
  if (!is.vector(Z_e) || !is.numeric(Z_e) || !all(Z_e %in% c(0, 1))) stop("Z_e must be a binary vector (0 or 1).")
  if (!is.vector(F_a) || !is.numeric(F_a) || !all(F_a %in% c(0, 1))) stop("F_a must be a binary vector (0 or 1).")
  if (!is.vector(ego_id_a) || !is.numeric(ego_id_a) || length(ego_id_a) != length(Y_a)) {
    stop("'ego_id_a' must be a numeric vector with the same length as Y_a.")
  }
  if (length(Y_e) != length(Z_e)) stop("Y_e and Z_e must have the same length (n_e).")
  if (length(Y_a) != length(F_a)) stop("Y_a and F_a must have the same length (n_a).")
  if (!is.numeric(kappa_vec)) stop("kappa_vec must be a numeric vector.")
  if (pz <= 0 | pz >= 1) stop("Invalid treatment probability 'pz' input.")

  n_e <- length(Y_e)

  # 1. Estimate outcome models using cross-fitting
  outcomes_res <- estimate_outcome_models_cf(
    Y_e = Y_e,
    Y_a = Y_a,
    X_e = X_e,
    X_a = X_a,
    Z_e = Z_e,
    F_a = F_a,
    ego_id_a = ego_id_a,
    reg_model_egos = reg_model_egos,
    reg_model_alters = reg_model_alters,
    formula_egos = formula_egos,
    formula_alters = formula_alters,
    n_folds = n_folds,
    ...
  )

  # 2. Bias-corrected IE estimates
  # 'IE_corrected' is data.table with columns: pi_param, ie_rd, ie_rd_var
  IE_corrected <- ie_aug_point_grid_(
    Y_a = Y_a,
    F_a = F_a,
    mu_a_1 = outcomes_res$mu_a_1,
    mu_a_0 = outcomes_res$mu_a_0,
    pi_list = pi_list_alter_ego,
    pz = pz,
    ego_id_a = ego_id_a,
    n_e = n_e,
    estimate_var = estimate_var
  )

  # 3. Bias-corrected DE estimates
  # 'DE_corrected' is data.table with columns: pi_param, kappa, de_rd, de_rd_var
  DE_corrected <- de_grid_multi_pi_kappa(
    Y_e = Y_e,
    Z_e = Z_e,
    mu_e_1 = outcomes_res$mu_e_1,
    mu_e_0 = outcomes_res$mu_e_0,
    pi_list = pi_list_ego_ego,
    kappa_vec = kappa_vec,
    pz = pz,
    estimate_var = estimate_var
  )

  return(list(
    IE_corrected = IE_corrected,
    DE_corrected = DE_corrected
  ))
}

