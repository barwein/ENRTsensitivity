

# Source + Packages --------------------------------------------------------------

# source("src/sensitivity_params.r")

# Indirect effects --------------------------------------------------------

#' @title Calculate augmented IE estimates for a grid of pi_params
#'
#' @param Y_a Vector of alter outcomes.
#' @param F_a Vector of alter observed exposures.
#' @param mu_a_1 Vector of cross-fit predictions E[Y|F=1,X].
#' @param mu_a_0 Vector of cross-fit predictions E[Y|F=0,X].
#' @param pi_list List of pi_vec sensitivity parameters.
#' @param pz Scalar, Pr(Z=1) among egos.
#' @param ego_id_a Vector mapping alters to ego indices (1 to n_e).
#' @param n_e Total number of egos.
#' @param folds_ids_a Vector indicating fold assignments for alters.
#'
#' @return A data.table with pi_param, ie_rd, ie_rd_var
#'
#' @keywords internal
#'
ie_aug_point_grid_ <- function(Y_a,
                               F_a,
                               mu_a_1 = NULL,
                               mu_a_0 = NULL,
                               pi_list,
                               pz,
                               ego_id_a = NULL,
                               n_e = NULL,
                               folds_ids_a = NULL) {
  # Input validation
  n_a <- length(Y_a)
  if (is.null(folds_ids_a)){
    folds_ids_a <- rep(1, n_a)
  }
  n_folds <- length(unique(folds_ids_a))

  if (is.null(mu_a_1) || is.null(mu_a_0)) {
    mu_a_1 <- rep(0, n_a)
    mu_a_0 <- rep(0, n_a)
    warning("mu_a_1 or mu_a_0 is NULL. Setting both to zero vectors.")
  } else{
    if (n_a != length(F_a) ||
        n_a != length(mu_a_1) ||
        n_a != length(mu_a_0)) {
      stop("Y_a, F_a, mu_a_1, and mu_a_0 must have the same length (n_a).")
    }
  }
  if (pz <= 0 | pz >= 1) {
    stop("Invalid treatment probability 'pz' input.")
  }

  results_list <- list()
  for (i in seq_along(pi_list)) {
    # pi_vec <- pi_list[[i]]
    pi_vec <- pi_list[[i]]$pi

    # Validate pi_vec
    if (any(pi_vec >= 1)) {
      warning("pi_vec contains values >= 1. IE is undefined.")
      results_list <- data.table::rbindlist(list(results_list,
                                                 data.table::data.table(
                                       pi_param = names(pi_list)[i],
                                       ie_rd = NA,
                                       ie_rd_var = NA
                                     )))
      next
    }

    # --- Fold-Specific Estimation ---
    est_k <- numeric(n_folds)
    var_k <- numeric(n_folds)
    n_k_vec <- numeric(n_folds) # Number of alters in fold k used for weighting


    for (k in unique(folds_ids_a)) {
      # Indices for this fold
      idx_a <- folds_ids_a == k

      # Subset data for this fold
      Y_a_k <- Y_a[idx_a]
      F_a_k <- F_a[idx_a]
      mu_a_1_k <- mu_a_1[idx_a]
      mu_a_0_k <- mu_a_0[idx_a]

      ego_id_a_k <- ego_id_a[idx_a]

      if(length(pi_vec) == 1) {
        pi_vec_k <- pi_vec
      } else{
        pi_vec_k <- pi_vec[idx_a]
      }

      # Fold size
      n_a_k <- sum(idx_a)
      n_k_vec[k] <- n_a_k

      # Determine unique egos in this fold to get n_e_k
      unique_egos_k <- unique(ego_id_a_k)

      if (n_a_k == 0) next

      # --- Point Estimate for Fold k ---
      # Bias adjustment terms
      term1 <- (Y_a_k - mu_a_1_k) * F_a_k / pz
      term2 <- (Y_a_k - mu_a_0_k) * (1 - F_a_k) / (1-pz)
      resid_diff <- mu_a_1_k - mu_a_0_k + term1 - term2
      weights_k <- (1 - pz) / (1 - pi_vec_k)
      weighted_resid_diff <- weights_k * resid_diff

      # Point estimates for fold k:
      est_k[k] <- mean(weighted_resid_diff, na.rm = TRUE) # This is IE_RD_k

      # Variance estimate for Fold k:
      # Sum over alters for each ego-network
      D_alter_i <- weights_k*(term1 - term2)
      T_i <- tapply(D_alter_i, ego_id_a_k, sum)
      mean_T_i <- mean(T_i, na.rm = TRUE)
      sum_sq_diff <- sum((T_i - mean_T_i) ^ 2, na.rm = TRUE)
      var_k[k] <- sum_sq_diff / n_a_k^2
    }

    # --- Aggregation ---
    # Global N
    N_total <- sum(n_k_vec)
    if (N_total != n_a){
      warning("Total number of alters across folds does not equal n_a.")
    }
    weights <- n_k_vec / N_total

    # Weighted Sum of Estimates
    ie_rd_agg <- sum(weights * est_k, na.rm = TRUE)
    ie_rd_var_agg <- sum((weights^2) * var_k, na.rm = TRUE)

    results_list <- data.table::rbindlist(list(results_list,
                                         data.table::data.table(
                                     pi_param = names(pi_list)[i],
                                     ie_rd = ie_rd_agg,
                                     ie_rd_var = ie_rd_var_agg
                                   )))
  }

  return(results_list)
}


# Direct effects ----------------------------------------------------------

#' @title Estimate bias-corrected DE estimates for a grid of (pi_param, kappa) combinations
#'
#' @param Y_e Vector of observed ego outcomes (n_e).
#' @param Z_e Vector of observed ego treatments (n_e).
#' @param mu_e_1 Vector of predictions E[Y|Z=1, X] (n_e).
#' @param mu_e_0 Vector of predictions E[Y|Z=0, X] (n_e).
#' @param pi_list A list of pi_vecs.
#' @param kappa_vec A vector of kappa values.
#' @param pz Scalar Pr(Z=1).
#' @param folds_ids_a Vector indicating fold assignments for alters.
#'
#' @return A data.table with columns: pi_param, kappa, de_rd, de_rd_var.
#'
#' @keywords internal
#'
de_grid_multi_pi_kappa <- function(Y_e,
                                   Z_e,
                                   mu_e_1 = NULL,
                                   mu_e_0 = NULL,
                                   pi_list,
                                   kappa_vec,
                                   pz,
                                   folds_ids_e = NULL) {

  # --- Input Checks ---
  n_e <- length(Y_e)
  if (is.null(folds_ids_e)){
    folds_ids_e <- rep(1, n_e)
  }
  n_folds <- length(unique(folds_ids_e))

  if(is.null(mu_e_1) || is.null(mu_e_0)) {
    warning("mu_e_1 or mu_e_0 is NULL. Setting both to zero vectors.")
    mu_e_1 <- rep(0, n_e)
    mu_e_0 <- rep(0, n_e)
  }  else{
    if (n_e != length(mu_e_1) ||
        n_e != length(mu_e_0)) {
      stop("Y_e, mu_e_1, and mu_e_0 must have the same length (n_e).")
    }
  }

  if (n_e != length(Z_e)) {
    stop("Y_e and Z_e must have the same length (n_e).")
  }
  # if (all(lapply(pi_list, function(p) { all(p >= 0) }) == FALSE)) {
  if (all(lapply(pi_list, function(p) { all(p$pi >= 0) }) == FALSE)) {
    stop("Some pi_vec contains negative values.")
  }

  # Two-dim grid of params
  grid_params <- expand.grid(pi_idx = seq_along(pi_list), kappa = kappa_vec)

  results_list <- list()

  for (r in 1:nrow(grid_params)){
    idx <- grid_params$pi_idx[r]
    k_val <- grid_params$kappa[r]
    # pi_vec <- pi_list[[idx]]
    pi_vec <- pi_list[[idx]]$pi
    rho_mat <- pi_list[[idx]]$rho

    est_k <- numeric(n_folds)
    var_k <- numeric(n_folds)
    n_k_vec <- numeric(n_folds)

    # Run over the folds
    for (k in 1:n_folds) {
      idx_e <- folds_ids_e == k
      n_e_k <- sum(idx_e)
      n_k_vec[k] <- n_e_k

      # Subset data
      Y_e_k <- Y_e[idx_e]
      Z_e_k <- Z_e[idx_e]
      mu_e_1_k <- mu_e_1[idx_e]
      mu_e_0_k <- mu_e_0[idx_e]
      if(length(pi_vec) == 1) {
        pi_vec_k <- pi_vec
      } else{
        pi_vec_k <- pi_vec[idx_e]
      }


      # --- Bias Adjustment (Ego level) ---
      term1 <- (Y_e_k - mu_e_1_k) * Z_e_k / pz
      term2 <- (Y_e_k - mu_e_0_k) * (1 - Z_e_k) / (1 - pz)
      resid_diff <- mu_e_1_k - mu_e_0_k + term1 - term2
      mean_pi_k <- mean(pi_vec_k, na.rm = TRUE)
      # weights_k <-  1 / (1 + pi_vec_k * (k_val - 1))
      weights_k <-  1 / (1 + mean_pi_k * (k_val - 1))
      # weighted_resid_diff <- weights_k * resid_diff

      # Point estimates for fold k:
      # est_k[k] <- mean(weighted_resid_diff, na.rm = TRUE) # This is DE_k
      est_k[k] <- weights_k*mean(resid_diff, na.rm = TRUE) # This is DE_k

      # Variance estimate for Fold k:
      # Sum over alters for each ego-network
      # D_ego_i <- weights_k*(term1 - term2)
      D_ego_i <- term1 - term2
      mean_D_i <- mean(D_ego_i, na.rm = TRUE)
      v_hat_k <- (D_ego_i - mean_D_i)^2
      # sum_sq_diff <- sum((D_ego_i - mean_D_i) ^ 2, na.rm = TRUE)
      sum_sq_diff <- sum(v_hat_k, na.rm = TRUE)
      var_neyman <- sum_sq_diff / (n_e_k^2)
      # var_k[k] <- sum_sq_diff / n_e_k^2

     # Compute VIF correction term using rho matrix
    if (!is.null(rho_mat) && is.matrix(rho_mat)) {

        # Subset rho for units in this fold
        # Note: We must subset both rows and cols to get the sub-graph of egos in this fold
        rho_k_sub <- rho_mat[idx_e, idx_e]

        # Calculate S matrix: Expected number of shared neighbors
        # S_ij = sum_k rho_ik * rho_jk
        S_mat <- rho_k_sub %*% t(rho_k_sub)

        # We only sum over i != j. Set diagonal to 0.
        # diag(S_mat) <- 0

        # Conservative bound term: S_ij * sqrt(v_i * v_j)
        sqrt_v <- sqrt(v_hat_k)

        # Calculate quadratic form: sqrt_v' * S * sqrt_v
        # This sums S_ij * sqrt(v_i) * sqrt(v_j) for all i, j
        cov_sum <- as.numeric(t(sqrt_v) %*% S_mat %*% sqrt_v)

        # Apply scaling factor: pz(1-pz) / N^2
        correction <- (pz * (1 - pz) / n_e_k^2) * cov_sum
      }

      # var_k[k] <- var_neyman + correction
      var_k[k] <- (var_neyman + correction)*(weights_k^2)
    }

    # Aggregation
    N_total <- sum(n_k_vec)
    weights <- n_k_vec / N_total

    de_rd_agg <- sum(weights * est_k, na.rm = TRUE)
    de_rd_var_agg <- sum((weights^2) * var_k, na.rm = TRUE)

    results_list <- data.table::rbindlist(list(
      results_list,
      data.table::data.table(
        pi_param = names(pi_list)[idx],
        kappa = k_val,
        de_rd = de_rd_agg,
        de_rd_var = de_rd_var_agg
      )
    ))
  }
  return(results_list)
}




