

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

  for (r in seq_len(nrow(grid_params))){
    idx <- grid_params$pi_idx[r]
    k_val <- grid_params$kappa[r]
    pi_vec <- pi_list[[idx]]$pi
    rho_mat <- pi_list[[idx]]$rho # n_e x n_e, diag = 0

    est_k <- numeric(n_folds)
    # var_k <- numeric(n_folds)
    var_neyman_k <- numeric(n_folds)
    u_k_vec <- numeric(n_folds)

    # Run over the folds
    for (k in 1:n_folds) {
      idx_e <- folds_ids_e == k
      n_e_k <- sum(idx_e)
      if (n_e_k == 0) next

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
      u_k <- n_e_k*(1 + mean_pi_k * (k_val - 1)) # effective sample size adjustment
      u_k_vec[k] <- u_k # store u_k value for later weighting of combined estimates

      # Point estimates for fold k:
      est_k[k] <- sum(resid_diff, na.rm = TRUE) / u_k # This is DE_k

      # Variance estimate for Fold k:
      # Sum over alters for each ego-network
      D_ego_i <- term1 - term2
      mean_D_i <- sum(D_ego_i, na.rm = TRUE) / u_k
      v_hat_k <- (D_ego_i - mean_D_i)^2
      var_neyman_k[k] <- sum(v_hat_k, na.rm = TRUE) / (u_k^2)
    }

      # # Contamination, within-fold, C_ij = min(1, rho_ij + xi_ij)  (A.10)
      # var_conta_k <- 0
      # if (!is.null(rho_mat) && is.matrix(rho_mat)) {
      #   rho_sub <- rho_mat[idx_e, idx_e, drop = FALSE]   # n_e_k x n_e_k
      #   xi_sub  <- rho_sub %*% t(rho_sub)                # xi_ij = sum_k rho_ik rho_jk
      #   # C_sub   <- pmin(rho_sub + xi_sub, 1)             # cap at 1
      #   C_sub   <- pmin(xi_sub, 1)             # cap at 1
      #   diag(C_sub) <- 0                                 # enforce i != j
      #   s_k       <- sqrt(v_hat_k)
      #   cov_sum_k <- as.numeric(t(s_k) %*% C_sub %*% s_k) # sum_{i!=j} C_ij sqrt(v_i v_j)
      #   var_conta_k <- cov_sum_k / (u_k^2)               # NO pz(1-pz)
      # }
      #
      # var_k[k] <- var_neyman_k + var_conta_k # var estimate in fold k

      # ------------------------------------------------------------------
      # Contamination-attributed variance (A.17), with C_ij (A.15) and
      # D_ij (A.14) computed within the current fold's egos only.
      #
      # V_Conta[q] = [ DE_adj[q] * (kappa - 1) ]^2 / u_e[q]^2
      #              * sum_{i != j in R_e ∩ S_q} ( C_ij + D_ij )
      #
      # Vanishes when kappa == 1 (covers the Naive case) or when rho == 0.
      # ------------------------------------------------------------------
    #   var_conta_k <- 0
    #   if (!is.null(rho_mat) && is.matrix(rho_mat)) {
    #
    #     # Full rho^e over all egos; enforce range and zero diagonal
    #     rho_full <- pmin(pmax(rho_mat, 0), 1)
    #     diag(rho_full) <- 0
    #
    #     # A_full[i,k] = 1 - pz * rho_ik, full n_e x n_e
    #     A_full   <- 1 - pz * rho_full
    #     logA_full <- log(A_full)                 # diag(logA_full) = 0
    #
    #     # log{ prod_{k != i,j} A_ik * A_jk } using factorization:
    #     #   = sum_k logA[i,k] + sum_k logA[j,k] - logA[i,i] - logA[j,j]
    #     #     - logA[i,j] - logA[j,i]
    #     #   = row_sum[i] + row_sum[j] - 2 * logA[i,j]   (logA symmetric, diag 0)
    #     log_row_sum_full <- rowSums(logA_full)
    #
    #     # Restrict to fold for (i, j) pairs
    #     log_row_sum_fold <- log_row_sum_full[idx_e]    # length n_e_k
    #     logA_fold        <- logA_full[idx_e, idx_e, drop = FALSE]
    #     rho_fold         <- rho_full[idx_e, idx_e, drop = FALSE]
    #     # Rows of rho restricted to fold egos but all columns (k spans R_e)
    #     rho_fold_rows    <- rho_full[idx_e, , drop = FALSE]   # n_e_k x n_e_tot
    #
    #     log_prod_a <- outer(log_row_sum_fold, log_row_sum_fold, "+") -
    #       2 * logA_fold
    #     prod_a <- exp(log_prod_a)              # n_e_k x n_e_k
    #
    #
    #     # log{ prod_{k in R_e, k != i,j} (1 - pz*(rho_ik + rho_jk + rho_ik*rho_jk)) }
    #     # Accumulate over ALL k in R_e via outer products on fold rows.
    #     log_first_full <- matrix(0, n_e_k, n_e_k)
    #     for (kk in seq_len(n_e_k)) {
    #       r_kk <- rho_fold_rows[, kk]          # rho_{i, kk} for i in fold
    #       B_kk <- 1 - pz * (outer(r_kk, r_kk, "+") + outer(r_kk, r_kk, "*"))
    #       # Guard against pathological negatives (shouldn't happen with rho in [0,1]
    #       # and pz in (0,1), but enforce numerically)
    #       B_kk[B_kk <= 0] <- .Machine$double.eps
    #       log_first_full <- log_first_full + log(B_kk)
    #     }
    #     # Remove the k = i and k = j terms.
    #     # When k = i: B_i[i,j] = 1 - p_z*(rho_ii + rho_ji + rho_ii*rho_ji)
    #     #                     = 1 - p_z*rho_ij = A_mat[i,j] (using rho_ii = 0).
    #     # Symmetric argument for k = j. So subtract 2 * logA[i,j].
    #     log_first  <- log_first_full - 2 * logA_fold
    #     first_prod <- exp(log_first)
    #
    #     # D_ij  (A.14), restricted to fold (i, j)
    #     D_mat <- pz^2 * rho_fold * (1 - rho_fold) * prod_a
    #     diag(D_mat) <- 0
    #
    #     # C_ij  (A.15), restricted to fold (i, j)
    #     coef_first  <- 1 - 2 * pz * rho_fold + (pz^2) * rho_fold
    #     coef_second <- (1 - pz * rho_fold)^2
    #     C_mat <- coef_first * first_prod - coef_second * prod_a
    #     diag(C_mat) <- 0
    #
    #
    #     sum_CD <- sum(C_mat) + sum(D_mat)                # diag already 0
    #     de_kappa_sq <- (est_k[k] * (k_val - 1))^2        # [DE_adj[q] * (kappa-1)]^2
    #     var_conta_k <- (de_kappa_sq / (u_k^2)) * sum_CD
    #   }
    #
    #   var_k[k] <- var_neyman_k + var_conta_k             # fold-k total variance
    #
    # }

    # --- Aggregate point estimate and Neyman variance across folds ---
    U_total       <- sum(u_k_vec)                  # u_e
    w             <- u_k_vec / U_total
    de_rd_agg     <- sum(w * est_k, na.rm = TRUE)
    var_neyman_agg <- sum((w^2) * var_neyman_k, na.rm = TRUE)

    # --- Global contamination-attributed variance (A.7) ---
    # Pairs (i, j) over ALL egos in R_e; products over k in R_e, k != i, j;
    # plug-in uses the aggregated estimate de_rd_agg.
    #
    # V_Conta = [ de_rd_agg * (kappa - 1) ]^2 / u_e^2
    #           * sum_{i != j in R_e} ( C_ij + D_ij ),
    # with C_ij and D_ij  as in the appendix.
    # Vanishes when kappa == 1 or rho == 0.
    var_conta_agg <- 0
    if (!is.null(rho_mat) && is.matrix(rho_mat) &&
        abs(k_val - 1) > 1e-12 && n_e >= 2) {

      # Full rho^e; enforce range and zero diagonal
      rho_full <- pmin(pmax(rho_mat, 0), 1)
      diag(rho_full) <- 0

      # A[i,k] = 1 - pz * rho_ik
      A_full    <- 1 - pz * rho_full
      logA_full <- log(A_full)                  # diag(logA_full) = 0

      # log{ prod_{k != i,j} A_ik * A_jk } via factorization
      log_row_sum <- rowSums(logA_full)
      # log_prod_a  <- outer(log_row_sum, log_row_sum, "+") - 2 * logA_full
      log_prod_a  <- outer(log_row_sum, log_row_sum, "+") - logA_full
      prod_a      <- exp(log_prod_a)            # n_e x n_e

      # log{ prod_{k != i,j} (1 - pz*(rho_ik + rho_jk + rho_ik*rho_jk)) }
      # Accumulate over all k in R_e.
      log_first_full <- matrix(0, n_e, n_e)
      for (kk in seq_len(n_e)) {
        r_kk <- rho_full[, kk]
        B_kk <- 1 - pz * (outer(r_kk, r_kk, "+") + outer(r_kk, r_kk, "*"))
        B_kk[B_kk <= 0] <- .Machine$double.eps  # numerical guard
        log_first_full <- log_first_full + log(B_kk)
      }
      # Remove k = i and k = j contributions (each equals log A_full[i,j],
      # since rho_ii = rho_jj = 0).
      log_first  <- log_first_full - logA_full
      # log_first  <- log_first_full - 2 * logA_full
      first_prod <- exp(log_first)

      # D_ij  (A.14)
      D_mat <- (pz^2) * rho_full * (1 - rho_full) * prod_a
      diag(D_mat) <- 0

      # C_ij  (A.15)
      coef_first  <- 1 - 2 * pz * rho_full + (pz^2) * rho_full
      coef_second <- (1 - pz * rho_full)^2
      C_mat <- coef_first * first_prod - coef_second * prod_a
      diag(C_mat) <- 0

      sum_CD      <- sum(C_mat) + sum(D_mat)
      de_kappa_sq <- (de_rd_agg * (k_val - 1))^2    # aggregate plug-in
      var_conta_agg <- (de_kappa_sq / (U_total^2)) * sum_CD
    }

    de_rd_var_agg <- var_neyman_agg + var_conta_agg

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




