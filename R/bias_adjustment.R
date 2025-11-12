

# Source + Packages --------------------------------------------------------------

# source("src/sensitivity_params.r")

# Indirect effects --------------------------------------------------------

#' @title Calculate augmented IE_RD and its cluster-robust variance for a single pi_vec
#'
#' @param Y_a Vector of alter outcomes.
#' @param F_a Vector of alter observed exposures.
#' @param mu_a_1 Vector of cross-fit predictions E[Y|F=1,X].
#' @param mu_a_0 Vector of cross-fit predictions E[Y|F=0,X].
#' @param pi_vec Vector of pi_i^a sensitivity parameters.
#' @param pz Scalar, Pr(Z=1) among egos.
#' @param ego_id_a Vector mapping alters to ego indices (1 to n_e).
#' @param n_e Total number of egos (for variance calculation).
#' @param estimate_var Logical. If TRUE, compute cluster-robust variance.
#'
#' @return A numeric vector c(ie_rd, ie_rd_var)
#'
#' @keywords internal
#'
ie_aug_point_one_param_ <- function(Y_a,
                                    F_a,
                                    mu_a_1,
                                    mu_a_0,
                                    pi_vec,
                                    pz,
                                    ego_id_a,
                                    n_e,
                                    estimate_var = FALSE) {

  n_a <- length(Y_a)

  # Input validation
  # Ensure pi_vec is valid
  if (any(pi_vec >= 1)) {
    warning("pi_vec contains values >= 1. IE estimate will be NA.")
    pi_vec[pi_vec >= 1] <- NA
  }

  # Weights for bias correction
  weights <- (1 - pz) / (1 - pi_vec)

  I_1 <- (F_a == 1)
  I_0 <- (F_a == 0)

  # Augmented estimator components
  term_1 <- (I_1 * (Y_a - mu_a_1)) / pz
  term_2 <- (I_0 * (Y_a - mu_a_0)) / (1 - pz)
  term_3 <- mu_a_1 - mu_a_0

  inside_sum <- term_1 - term_2 + term_3

  # Element-wise product: alpha_i^a * phi_i^a
  alpha_phi_a <- weights * inside_sum
  # Final estimator for RD
  ie_rd_ <- mean(alpha_phi_a, na.rm = TRUE)

  # --- Variance Estimation ---
  # Estimate variance with cluster-robust empirical variance estimator
  ie_rd_var <- NA_real_

  if (estimate_var) {
    if (is.null(ego_id_a) || is.null(n_e)) {
      stop("ego_id_a and n_e must be provided to estimate cluster-robust variance for IE.")
    }
    if (n_e <= 1) {
      warning("Cannot compute variance with n_e <= 1. Setting var to NA.")
    } else {
      # Create T_i for all egos (i in R_e)
      # Egos with no alters will have T_i = 0
      T_full <- numeric(n_e)

      # Sum (alpha_j^a * phi_j^a) for alters grouped by their ego
      # This creates a named vector where names are ego indices
      T_i_vec <- tapply(alpha_phi_a, ego_id_a, sum, na.rm = TRUE)

      # Get the ego indices that actually have alters
      ego_indices_with_alters <- as.numeric(names(T_i_vec))

      # Place the sums into the full vector
      T_full[ego_indices_with_alters] <- T_i_vec

      # T_bar = (1/n_e) * sum(T_i)
      T_bar <- mean(T_full)

      # Sum of squared differences
      sum_sq_diff <- sum((T_full - T_bar) ^ 2, na.rm = TRUE)

      # Variance estimator from appendix
      ie_rd_var <- (n_e / (n_a * n_a * (n_e - 1))) * sum_sq_diff
    }
  }

  return(c(ie_rd = ie_rd_, ie_rd_var = ie_rd_var))
}


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
#' @param estimate_var Logical. If TRUE, compute cluster-robust variance.
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
                               estimate_var = FALSE) {
  # Input validation
  n_a <- length(Y_a)
  if (is.null(mu_a_1) || is.null(mu_a_0)) {
    mu_a_1 <- rep(0, n_a)
    mu_a_0 <- rep(0, n_a)
    warning("mu_a_1 or mu_a_0 is NULL. Setting both to zero vectors.")
  }
  else{
    if (n_a != length(F_a) ||
        n_a != length(mu_a_1) ||
        n_a != length(mu_a_0)) {
      stop("Y_a, F_a, mu_a_1, and mu_a_0 must have the same length (n_a).")
    }
  }
  if (pz <= 0 | pz >= 1) {
    stop("Invalid treatment probability 'pz' input.")
  }

  ie_grid_mat <- vapply(pi_list, function(pi_v) {
    # Ensure pi_v has correct length if it's not scalar
    if (length(pi_v) > 1 && length(pi_v) != n_a) {
      stop(paste("Heterogeneous pi_vec length", length(pi_v),
                 "does not match alter sample size", n_a))
    }

    out_vec <- ie_aug_point_one_param_(
      Y_a = Y_a,
      F_a = F_a,
      mu_a_1 = mu_a_1,
      mu_a_0 = mu_a_0,
      pi_vec = pi_v,
      pz = pz,
      ego_id_a = ego_id_a,
      n_e = n_e,
      estimate_var = estimate_var
    )
    out_vec
  }, FUN.VALUE = numeric(2)) # Expecting a vector of length 2 (ie_rd, ie_rd_var)

  res_dt <- data.table(
    pi_param = as.numeric(names(pi_list)),
    ie_rd = ie_grid_mat[1,],      # First row is ie_rd
    ie_rd_var = ie_grid_mat[2,] # Second row is ie_rd_var
  )
  return(res_dt)
}


# Direct effects ----------------------------------------------------------


#' @title Calculate bias-corrected DE_RD for a single (pi, kappa) combination
#'
#' @description
#' Implements the augmented randomization-based estimator:
#' DE_aug = (1/n_e) * sum[ w_i * phi_i ]
#' where w_i = 1 / (1 + pi_i * (kappa - 1))
#' and   phi_i = [ I(Z=1)(Y-m1)/pz - I(Z=0)(Y-m0)/(1-pz) + (m1 - m0) ]
#'
#' @param Y_e Vector of observed ego outcomes.
#' @param Z_e Vector of observed ego treatments.
#' @param mu_e_1 Vector of predictions m^e_i(1) = E[Y|Z=1, X].
#' @param mu_e_0 Vector of predictions m^e_i(0) = E[Y|Z=0, X].
#' @param pi_vec Vector of pi_i^e values.
#' @param kappa_ Scalar kappa value.
#' @param pz Scalar Pr(Z=1).
#' @param estimate_var Logical. If TRUE, compute empirical variance.

#' @return A numeric vector `c(de_rd, de_rd_var)`
#'
#' @keywords internal
#'
de_point_one_pi_kappa <- function(Y_e,
                                  Z_e,
                                  mu_e_1,
                                  mu_e_0,
                                  pi_vec,
                                  kappa_,
                                  pz,
                                  estimate_var = FALSE) {

  # --- 1. Calculate the augmented difference term (A_i) ---
  I_1 <- (Z_e == 1)
  I_0 <- (Z_e == 0)

  term_1 <- (I_1 * (Y_e - mu_e_1)) / pz
  term_2 <- (I_0 * (Y_e - mu_e_0)) / (1 - pz)
  term_3 <- mu_e_1 - mu_e_0

  # Vector of augmented differences, one per ego
  aug_diff <- term_1 - term_2 + term_3

  # --- 2. Calculate the DE weight ---
  # Ensure pi_vec is broadcastable if it's a scalar
  if (length(pi_vec) == 1) {
    pi_vec <- rep(pi_vec, length(Y_e))
  }

  de_weight <- 1 / (1 + pi_vec * (kappa_ - 1))

  if (any(is.infinite(de_weight)) || any(is.na(de_weight))) {
    warning("DE weights are non-finite. Check pi_vec and kappa values. Setting weights to zero.")
    de_weight[is.infinite(de_weight) | is.na(de_weight)] <- 0
  }

  # --- 3. Calculate final estimator (alpha_i^e * phi_i^e) ---
  alpha_phi_e <- de_weight * aug_diff
  de_rd_ <- mean(alpha_phi_e, na.rm = TRUE)

  # --- 4. Variance Estimation ---
  de_rd_var <- NA_real_
  n_e <- length(Y_e)
  if (estimate_var) {
    if (n_e <= 1) {
      warning("Cannot compute variance with n_e <= 1. Setting var to NA.")
    } else {
      # alpha_phi_e - DE_aug
      sum_sq_diff <- sum((alpha_phi_e - de_rd_) ^ 2, na.rm = TRUE)

      # Variance estimator from appendix
      de_rd_var <- (1 / (n_e * (n_e - 1))) * sum_sq_diff
    }
  }

  return(c(de_rd = de_rd_, de_rd_var = de_rd_var))
}

#' @title Estimate bias-corrected DE for a single pi_param and a grid of kappa values
#'
#' @param Y_e Vector of observed ego outcomes.
#' @param Z_e Vector of observed ego treatments.
#' @param mu_e_1 Vector of predictions E[Y|Z=1, X].
#' @param mu_e_0 Vector of predictions E[Y|Z=0, X].
#' @param pi_vec Vector of pi_i^e values.
#' @param kappa_vec_ Vector of kappa values to iterate over.
#' @param pz Scalar Pr(Z=1).
#' @param estimate_var Logical. If TRUE, compute empirical variance.
#'
#' @return A data.table with columns: kappa, de_rd, de_rd_var
#'
#' @keywords internal
#'
de_grid_one_pi_multi_kappa <- function(Y_e,
                                       Z_e,
                                       mu_e_1,
                                       mu_e_0,
                                       pi_vec,
                                       kappa_vec_,
                                       pz,
                                       estimate_var = FALSE) {

  # Ensure pi_vec is a numeric vector
  pi_vec_numeric <- unlist(pi_vec)

  # Create a named list of kappa values for vapply
  kappa_list <- as.list(kappa_vec_)
  names(kappa_list) <- round(kappa_vec_, 3)

  # Iterate over all kappa values, expecting 2 return values (mean, var)
  de_kappa_mat <- vapply(kappa_list, function(kappa_) {
    out_vec <- de_point_one_pi_kappa(
      Y_e = Y_e,
      Z_e = Z_e,
      mu_e_1 = mu_e_1,
      mu_e_0 = mu_e_0,
      pi_vec = pi_vec_numeric,
      kappa_ = kappa_,
      pz = pz,
      estimate_var = estimate_var
    )
    out_vec
  },
  FUN.VALUE = numeric(2)) # Expecting vector of length 2

  res_dt <- data.table(
    kappa = names(kappa_list),
    de_rd = de_kappa_mat[1,],      # First row is de_rd
    de_rd_var = de_kappa_mat[2,]  # Second row is de_rd_var
  )
  return(res_dt)
}


#' @title Estimate bias-corrected DE estimates for a grid of (pi_param, kappa) combinations
#'
#' @param Y_e Vector of observed ego outcomes (n_e).
#' @param Z_e Vector of observed ego treatments (n_e).
#' @param mu_e_1 Vector of predictions E[Y|Z=1, X] (n_e).
#' @param mu_e_0 Vector of predictions E[Y|Z=0, X] (n_e).
#' @param pi_list A list of pi_vecs.
#' @param kappa_vec A vector of kappa values.
#' @param pz Scalar Pr(Z=1).
#' @param estimate_var Logical. If TRUE, compute empirical variance.
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
                                   estimate_var = FALSE) {

  # --- Input Checks ---
  n_e <- length(Y_e)

  if(is.null(mu_e_1) || is.null(mu_e_0)) {
    warning("mu_e_1 or mu_e_0 is NULL. Setting both to zero vectors.")
    mu_e_1 <- rep(0, n_e)
    mu_e_0 <- rep(0, n_e)
  }
  else{
    if (n_e != length(mu_e_1) ||
        n_e != length(mu_e_0)) {
      stop("Y_e, mu_e_1, and mu_e_0 must have the same length (n_e).")
    }
  }

  if (n_e != length(Z_e)) {
    stop("Y_e and Z_e must have the same length (n_e).")
  }
  if (all(lapply(pi_list, function(p) { all(p >= 0) }) == FALSE)) {
    stop("Some pi_vec contains negative values.")
  }

  # Iterate over all pi_list entries
  de_list <- lapply(pi_list, function(pi_v) {
    de_grid_one_pi_multi_kappa(
      Y_e = Y_e,
      Z_e = Z_e,
      mu_e_1 = mu_e_1,
      mu_e_0 = mu_e_0,
      pi_vec = pi_v,
      kappa_vec_ = kappa_vec,
      pz = pz,
      estimate_var = estimate_var
    )
  })

  res_dt <- rbindlist(de_list, idcol = "pi_param")
  res_dt$pi_param <- as.numeric(res_dt$pi_param)
  res_dt$kappa <- as.numeric(res_dt$kappa)

  return(res_dt)
}




