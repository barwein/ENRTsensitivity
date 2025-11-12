
# Load required source files
# source("src/bias_adjustment.r")
# source("src/outcomes_models.R")
# source("src/sa_single_iter.R")
#
# # Load required libraries
# library(data.table)
# library(parallel)

#' @title Wrapper for Sensitivity Analysis with Non-Parametric Bootstrap
#'
#' @description
#' This function wraps the single-iteration sensitivity analysis (`SA_one_iter`)
#' in a non-parametric bootstrap procedure to estimate the variance of the
#' bias-corrected indirect and direct effects.
#' The bootstrap procedure resamples ego-networks (i.e., egos and all their
#' recruited alters) with replacement. The variance is estimated by computing
#' the sample variance of the point estimates obtained from the bootstrap samples.
#'
#' @param Y_e A numeric vector of outcomes for the egos.
#' @param Y_a A numeric vector of outcomes for the alters.
#' @param X_e A numeric matrix of covariates for the egos.
#' @param X_a A numeric matrix of covariates for the alters.
#' @param Z_e A binary numeric vector (0 or 1) of treatment assignments for egos.
#' @param F_a A binary numeric vector (0 or 1) of `observed` exposures for alters.
#' @param ego_id_a A numeric vector mapping each alter to their ego's index
#'        (an integer from 1 to n_e).
#' @param reg_model_egos A function for the egos' outcome regression model, e.g., `glm`, `lm`.
#' @param reg_model_alters A function for the alters' outcome regression model, e.g., `glm`, `lm`.
#' @param formula_egos A formula for the egos' regression model.
#' @param formula_alters A formula for the alters' regression model.
#' @param pi_list_ego_ego A list of exposure probabilities (`pi_i^e`) for egos.
#' @param pi_list_alter_ego A list of exposure probabilities (`pi_i^a`) for alters.
#' @param kappa_vec A numeric vector of values for the kappa sensitivity parameter.
#' @param pz The known probability of an ego being assigned to treatment, Pr(Z=1).
#' @param B An integer specifying the number of bootstrap iterations.
#' @param n_folds The number of folds to use for cross-fitting outcome models
#'        in `SA_one_iter`.
#' @param n_cores An integer specifying the number of CPU cores to use for
#'        parallel processing.
#' @param verbose A logical indicating whether to print progress messages.
#' @param ... Additional arguments passed to the regression model functions.
#'
#' @return A list containing two data.tables:
#' \describe{
#'   \item{IE_results}{A data.table with point estimates (`ie_rd`),
#'         empirical variance (`ie_rd_var`), and bootstrap variance
#'         (`ie_rd_boot_var`) for each sensitivity parameter.}
#'   \item{DE_results}{A data.table with point estimates (`de_rd`),
#'         empirical variance (`de_rd_var`), and bootstrap variance
#'         (`de_rd_boot_var`) for each combination of sensitivity parameters.}
#' }
#' @keywords internal
#'
run_sensitivity_bootstrap <- function(Y_e,
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
                                      pz = 0.5,
                                      B = 500,
                                      n_folds = 2,
                                      n_cores = 1,
                                      verbose = TRUE,
                                      ...) {

  # --- Input Validation ---
  if (!is.vector(ego_id_a) || !is.numeric(ego_id_a) || length(ego_id_a) != length(Y_a)) {
    stop("'ego_id_a' must be a numeric vector with the same length as the number of alters.")
  }
  if (max(ego_id_a) > length(Y_e)) {
    stop("An ID in 'ego_id_a' references an ego index that does not exist.")
  }

  n_e <- length(Y_e)
  n_a <- length(Y_a)

  # --- 1. Point Estimates & Empirical Variance (Full Sample) ---

  if(verbose) message("Calculating point estimates and empirical variance on the full dataset...")

  full_sample_results <- SA_one_iter(
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
    pi_list_ego_ego = pi_list_ego_ego,
    pi_list_alter_ego = pi_list_alter_ego,
    kappa_vec = kappa_vec,
    pz = pz,
    n_folds = n_folds,
    estimate_var = TRUE, # Get empirical variance
    ...
  )

  IE_point_estimates <- full_sample_results$IE_corrected
  DE_point_estimates <- full_sample_results$DE_corrected

  # Ensure kappa_vec is defined from the full run
  kappa_vec = unique(DE_point_estimates$kappa)

  # --- 2. Non-Parametric Bootstrap ---

  # Pre-split alter indices by ego for efficient resampling
  alter_indices_by_ego <- split(1:n_a, ego_id_a)

  # Define a helper function to perform a single bootstrap iteration.
  bootstrap_iteration <- function(iter) {
    # Resample ego indices with replacement
    boot_ego_indices <- sample(1:n_e, size = n_e, replace = TRUE)

    # Get the list of alter indices corresponding to the sampled egos
    boot_alter_list <- alter_indices_by_ego[boot_ego_indices]

    # Get the vector of alter indices
    boot_alter_indices <- unlist(boot_alter_list, use.names = FALSE)

    # Create the new ego_id_a vector, mapping alters to their new
    # *bootstrap* ego index (from 1 to n_e)
    n_alters_per_boot_ego <- sapply(boot_alter_list, length)
    ego_id_a_boot <- rep(1:n_e, times = n_alters_per_boot_ego)

    # Create bootstrap datasets
    Y_e_boot <- Y_e[boot_ego_indices]
    Z_e_boot <- Z_e[boot_ego_indices]
    X_e_boot <- if (!is.null(X_e)) X_e[boot_ego_indices, , drop = FALSE] else NULL

    Y_a_boot <- Y_a[boot_alter_indices]
    F_a_boot <- F_a[boot_alter_indices]
    X_a_boot <- if (!is.null(X_a)) X_a[boot_alter_indices, , drop = FALSE] else NULL

    # Subset heterogeneous pi vectors if they are not scalars
    pi_list_ego_boot <- lapply(pi_list_ego_ego, function(pi_vec) {
      if (length(pi_vec) > 1) pi_vec[boot_ego_indices] else pi_vec
    })

    pi_list_alter_boot <- lapply(pi_list_alter_ego, function(pi_vec) {
      if (length(pi_vec) > 1) pi_vec[boot_alter_indices] else pi_vec
    })

    # Run SA for the bootstrap sample, return NULL on error
    tryCatch({
      SA_one_iter(
        Y_e = Y_e_boot,
        Y_a = Y_a_boot,
        X_e = X_e_boot,
        X_a = X_a_boot,
        Z_e = Z_e_boot,
        F_a = F_a_boot,
        ego_id_a = ego_id_a_boot, # Use the re-indexed ego IDs
        reg_model_egos = reg_model_egos,
        reg_model_alters = reg_model_alters,
        formula_egos = formula_egos,
        formula_alters = formula_alters,
        pi_list_ego_ego = pi_list_ego_boot,
        pi_list_alter_ego = pi_list_alter_boot,
        kappa_vec = kappa_vec,
        pz = pz,
        n_folds = n_folds,
        estimate_var = FALSE, # Do not need variance from the bootstrap sample
        ...
      )
    }, error = function(e) {
      warning(paste("Bootstrap iteration failed:", e$message))
      return(NULL)
    })
  }

  if(verbose) message(paste0("Starting ", B, " bootstrap iterations on ", n_cores, " core(s)..."))

  # Choose the execution method based on n_cores
  if (n_cores > 1 && .Platform$OS.type != "windows") {
    boot_results_list <- mclapply(1:B, bootstrap_iteration, mc.cores = n_cores)
  } else {
    if (n_cores > 1 && .Platform$OS.type == "windows") {
      warning("mclapply is not supported on Windows. Using sequential 'lapply'.")
    }
    boot_results_list <- lapply(1:B, bootstrap_iteration)
  }

  # --- 3. Aggregate Bootstrap Results ---

  if(verbose) message("Aggregating bootstrap results...")

  # Filter out any NULL results from failed iterations
  boot_results_list <- boot_results_list[!sapply(boot_results_list, is.null)]

  # Separate the IE and DE results into their own lists
  boot_ie_dt <- rbindlist(lapply(boot_results_list, `[[`, "IE_corrected"))
  boot_de_dt <- rbindlist(lapply(boot_results_list, `[[`, "DE_corrected"))

  # Summarize IE results: calculate variance of the point estimates
  IE_summary <- boot_ie_dt[, .(
    ie_rd_boot_var = var(ie_rd, na.rm = TRUE)
  ), by = pi_param]

  # Summarize DE results: calculate variance of the point estimates
  DE_summary <- boot_de_dt[, .(
    de_rd_boot_var = var(de_rd, na.rm = TRUE)
  ), by = .(pi_param, kappa)]

  # --- 4. Combine results ---

  # Combine IE (merges full sample results with bootstrap variance)
  IE_results <- merge(IE_point_estimates, IE_summary, by = "pi_param")

  # Combine DE (merges full sample results with bootstrap variance)
  DE_results <- merge(DE_point_estimates, DE_summary, by = c("pi_param", "kappa"))


  # --- 5. Return Final Results ---

  return(list(
    IE_results = IE_results,
    DE_results = DE_results
  ))
}

