# Source dependencies
# source("src/bias_adjustment.r")
# source("src/outcomes_models.R")
# source("src/sa_single_iter.R")
# source("src/sensitivity_params.r")
#
# library(data.table)
# library(parallel)

#' @title Perform Probabilistic Bias Analysis (PBA) for contamination in ENRT
#' @description
#' This function performs a PBA as described in Section 3.4 of the paper.
#' It integrates over uncertainty in sensitivity
#' parameters and sampling uncertainty. It uses a single Monte Carlo loop to
#' efficiently generate distributions for two types of uncertainty:
#' 1.  **Bias-Only Uncertainty:** Reflects uncertainty from the sensitivity
#'     parameters, holding the dataset fixed.
#' 2.  **Total Uncertainty:** Reflects both sampling and sensitivity
#'     parameter uncertainty.
#'
#'
#' PBA uses a Normal approximation method. In each
#'     iteration, it draws a sample from the sensitivity parameter priors,
#'     calculates the bias-corrected point estimate and its empirical variance
#'     on the *full dataset*, and then draws from the resulting normal
#'     approximation `N(estimate, variance)` to represent total uncertainty.
#'
#' @param Y_e A numeric vector of outcomes for the egos.
#' @param Y_a A numeric vector of outcomes for the alters.
#' @param X_e A numeric matrix of covariates for the egos.
#' @param X_a A numeric matrix of covariates for the alters.
#' @param Z_e A binary numeric vector (0 or 1) of treatment assignments for egos.
#' @param F_a A binary numeric vector (0 or 1) of `observed` exposures for alters.
#' @param ego_id_a A numeric vector mapping each alter to their ego's index.
#'        Indices should correspond to the rows in the ego datasets (Y_e, X_e).
#' @param augmented Logical whether to use augmented estimators
#' @param reg_model_egos A function for the egos' outcome regression model, e.g., `glm`.
#' @param reg_model_alters A function for the alters' outcome regression model, e.g., `glm`.
#' @param formula_egos A formula for the egos' regression model.
#' @param formula_alters A formula for the alters' regression model.
#' @param B An integer specifying the number of Monte Carlo iterations.
#' @param pz The known probability of an ego being assigned to treatment, Pr(Z=1).
#' @param probs A numeric vector (length=`2`) of probabilities for quantiles.
#' @param n_cores Number of cores for parallel execution.
#' @param n_folds An integer specifying the number of folds for cross-fitting.
#' @param verbose A logical indicating whether to print progress messages.
#' @param prior_func_ie A function that returns a single sampled value for the IE
#'        sensitivity parameter (e.g., one 'rho_a' or 'm_a').
#' @param pi_func_ie The R function to calculate the IE pi vector (e.g., `pi_homo`).
#' @param pi_args_ie A list of all arguments for `pi_func_ie`.
#' @param pi_param_name_ie The string name of the argument in `pi_func_ie`.
#' @param prior_func_de A function that returns a list with two named elements:
#'        'pi_param' (sampled DE pi parameter) and 'kappa'.
#' @param pi_func_de The R function to calculate the DE pi vector (e.g., `pi_homo`).
#' @param pi_args_de A list of arguments for `pi_func_de`.
#' @param pi_param_name_de The string name of the pi parameter in `pi_func_de`.
#'
#' @param ... Additional arguments passed to the regression model functions.
#'
#' @return A list containing two data.tables: `IE_results` and `DE_results`.
#'         Each table provides the mean and quantiles for both the
#'         'bias_only' and 'total' uncertainty distributions.
#'
#' @importFrom data.table := .N .SD data.table rbindlist
#' @importFrom stats qnorm rnorm quantile predict
#' @importFrom parallel mclapply
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # --- 1. Simulate HPTN-037-like data (same as enrt_sa) ---
#' set.seed(42)
#' n_e <- 150
#' n_a <- 263
#' pz <- 0.5
#'
#' # (Egos)
#' X_e <- matrix(rnorm(n_e * 2), n_e, 2, dimnames = list(NULL, c("X1", "X2")))
#' Z_e <- rbinom(n_e, 1, pz)
#' Y_e <- rbinom(n_e, 1, plogis(0.3 + 0.1 * Z_e - 0.2 * X_e[,1]))
#'
#' # (Alters)
#' ego_id_a <- sample(1:n_e, n_a, replace = TRUE)
#' X_a <- matrix(rnorm(n_a * 2), n_a, 2, dimnames = list(NULL, c("X1", "X2")))
#' F_a <- Z_e[ego_id_a] # Observed exposure
#' Y_a <- rbinom(n_a, 1, plogis(0.4 - 0.1 * F_a + 0.1 * X_a[,2]))
#'
#' # --- 2. Define Prior Distributions and Arguments ---
#'
#' # Priors for IE (Homogeneous m_a, Example 2)
#' # Use a Poisson(200) prior for m_a
#' prior_ie_func <- function() {
#'   round(rpois(1, 200)) # Sample a single m_a value
#' }
#'
#' # Arguments for pi_homo (for alters)
#' pi_args_ie_list <- list(
#'   n_e = n_e,
#'   n_a = n_a,
#'   type = "alter",
#'   pz = pz
#' )
#'
#' # Priors for DE (Homogeneous m_e and kappa)
#' # Use Poisson(150) for m_e and Uniform(1, 3) for kappa.
#' prior_de_func <- function() {
#'   list(
#'     pi_param = round(rpois(1, 150)), # Sampled m_e
#'     kappa = runif(1, 1, 3)           # Sampled kappa
#'   )
#' }
#'
#' # Arguments for pi_homo (for egos)
#' pi_args_de_list <- list(
#'   n_e = n_e,
#'   type = "ego",
#'   pz = pz
#' )
#'
#' # --- 3. Run the Probabilistic Bias Analysis ---
#'
#' # Using augmented estimators (logistic regression)
#'
#' pba_results <- enrt_pba(
#'   Y_e = Y_e, Y_a = Y_a, X_e = X_e, X_a = X_a,
#'   Z_e = Z_e, F_a = F_a, ego_id_a = ego_id_a,
#'   augmented = TRUE,
#'   reg_model_egos = glm, reg_model_alters = glm,
#'   formula_egos = Y ~ Z + X1 + X2,
#'   formula_alters = Y ~ F + X1 + X2,
#'   B = 1000, # Recommend 1e3 to 1e4 for full analysis
#'   pz = pz,
#'   n_cores = 1, # modify
#'   n_folds = 2,
#'   prior_func_ie = prior_ie_func,
#'   pi_func_ie = pi_homo,
#'   pi_args_ie = pi_args_ie_list,
#'   pi_param_name_ie = "m_vec", # Name of arg to pass sampled value to
#'   prior_func_de = prior_de_func,
#'   pi_func_de = pi_homo,
#'   pi_args_de = pi_args_de_list,
#'   pi_param_name_de = "m_vec", # Name of arg to pass pi_param to
#'   family = binomial(link = "logit") # Passed to glm()
#' )
#'
#' # View the 95% intervals for bias-only and total uncertainty
#' print(pba_results$IE_results)
#' print(pba_results$DE_results)
#' }
#'
enrt_pba <- function(Y_e,
                     Y_a,
                     X_e = NULL,
                     X_a = NULL,
                     Z_e,
                     F_a,
                     ego_id_a,
                     augmented = FALSE,
                     reg_model_egos = NULL,
                     reg_model_alters= NULL,
                     formula_egos = NULL,
                     formula_alters = NULL,
                     # PBA params
                     B = 1e3,
                     pz = 0.5,
                     probs = c(0.025, 0.975),
                     n_cores = 1,
                     n_folds = 2,
                     verbose = TRUE,
                     # IE Prior / Pi specifications
                     prior_func_ie,
                     pi_func_ie,
                     pi_args_ie,
                     pi_param_name_ie,
                     # DE Prior / Pi specifications
                     prior_func_de,
                     pi_func_de,
                     pi_args_de,
                     pi_param_name_de,
                     ...) {

  # --- Input Validation ---
  n_e <- length(Y_e)
  n_a <- length(Y_a)

  if (!is.function(prior_func_ie)) stop("'prior_func_ie' must be a function.")
  if (!is.function(pi_func_ie)) stop("'pi_func_ie' must be a function.")
  if (!is.list(pi_args_ie)) stop("'pi_args_ie' must be a list.")
  if (!is.character(pi_param_name_ie) || length(pi_param_name_ie) != 1) {
    stop("'pi_param_name_ie' must be a single string.")
  }

  if (!is.function(prior_func_de)) stop("'prior_func_de' must be a function.")
  if (!is.function(pi_func_de)) stop("'pi_func_de' must be a function.")
  if (!is.list(pi_args_de)) stop("'pi_args_de' must be a list.")
  if (!is.character(pi_param_name_de) || length(pi_param_name_de) != 1) {
    stop("'pi_param_name_de' must be a single string.")
  }


  # --- Helper: Summarization Function (RD Only) ---
  summarize_pba_results <- function(results, type) {
    # Extract the 'IE_corrected' or 'DE_corrected' data.table from each list element
    ie_dt <- data.table::rbindlist(lapply(results, `[[`, "IE_corrected"))
    de_dt <- data.table::rbindlist(lapply(results, `[[`, "DE_corrected"))

    ie_summary <- ie_dt[, .(
      ie_rd_mean = mean(ie_rd, na.rm = TRUE),
      ie_rd_q_low = quantile(ie_rd, probs = probs[1], na.rm = TRUE),
      ie_rd_q_high = quantile(ie_rd, probs = probs[2], na.rm = TRUE)
    )]

    de_summary <- de_dt[, .(
      de_rd_mean = mean(de_rd, na.rm = TRUE),
      de_rd_q_low = quantile(de_rd, probs = probs[1], na.rm = TRUE),
      de_rd_q_high = quantile(de_rd, probs = probs[2], na.rm = TRUE)
    )]

    return(list(
      IE_results = data.table::data.table(uncertainty_type = type, ie_summary),
      DE_results = data.table::data.table(uncertainty_type = type, de_summary)
    ))
  }


  # --- Estimate the outcome models once ---
  # Estimate outcome models using cross-fitting
  if (augmented){
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
      # family = binomial(link = "logit") # Passed to glm()
    )
  } else{ # Non-augmented case (no cross-fitting)
    outcomes_res <- list(
      mu_a_1 = rep(0, n_a), mu_a_0 = rep(0, n_a),
      mu_e_1 = rep(0, n_e), mu_e_0 = rep(0, n_e),
      folds_ids_e = rep(1, n_e), # Single fold containing all indices
      folds_ids_a = rep(1, n_a)
    )
  }



  # --- Iteration Function  ---
  iteration_function <- function(iter) {
    # 1. Sample sensitivity parameters *once* for this iteration
    s_a <- prior_func_ie()
    s_de <- prior_func_de()

    # --- Calculate Bias-Only Estimates & Variance (Full Data D, Sampled s) ---
    res_full_sample <- tryCatch({
      # IE Pi List (full data)
      pi_args_alter_full <- pi_args_ie
      pi_args_alter_full[[pi_param_name_ie]] <- s_a
      pi_list_alter_full <- do.call(pi_func_ie, pi_args_alter_full)

      # DE Pi List (full data)
      pi_args_ego_full <- pi_args_de
      pi_args_ego_full[[pi_param_name_de]] <- s_de$pi_param
      pi_list_ego_full <- do.call(pi_func_de, pi_args_ego_full)

      SA_one_iter(
        Y_e = Y_e, Y_a = Y_a,
        X_e = X_e, X_a = X_a,
        Z_e = Z_e, F_a = F_a,
        ego_id_a = ego_id_a,
        mu_a_1 = outcomes_res$mu_a_1,
        mu_a_0 = outcomes_res$mu_a_0,
        mu_e_1 = outcomes_res$mu_e_1,
        mu_e_0 = outcomes_res$mu_e_0,
        folds_ids_a = outcomes_res$folds_ids_a,
        folds_ids_e = outcomes_res$folds_ids_e,
        pi_list_ego_ego = pi_list_ego_full,
        pi_list_alter_ego = pi_list_alter_full,
        kappa_vec = s_de$kappa,
        pz = pz,
        n_folds = n_folds,
        ...
      )
    }, error = function(e) {
      warning(paste("Full-sample calculation failed in iter", iter, ":", e$message))
      return(NULL)
    })

    if (is.null(res_full_sample)) {
      return(list(bias_only = NULL, total = NULL))
    }

    # --- A. 'bias_only' results are just the point estimates from the full sample ---
    res_bias_only <- res_full_sample

    # --- B. 'total' results are sampled from N(mean, var) ---

    # Sample for IE
    ie_res <- res_full_sample$IE_corrected
    ie_rd_total <- NA_real_
    # Check for valid mean and non-negative variance before sampling
    if (!is.na(ie_res$ie_rd) && !is.na(ie_res$ie_rd_var) && ie_res$ie_rd_var >= 0) {
      ie_rd_total <- rnorm(1, mean = ie_res$ie_rd, sd = sqrt(ie_res$ie_rd_var))
    }
    # Create the data.table for this iteration's *total* estimate
    res_total_ie <- data.table::data.table(
      pi_param = ie_res$pi_param,
      ie_rd = ie_rd_total
    )

    # Sample for DE
    de_res <- res_full_sample$DE_corrected
    de_rd_total <- NA_real_
    if (!is.na(de_res$de_rd) && !is.na(de_res$de_rd_var) && de_res$de_rd_var >= 0) {
      de_rd_total <- rnorm(1, mean = de_res$de_rd, sd = sqrt(de_res$de_rd_var))
    }
    res_total_de <- data.table(
      pi_param = de_res$pi_param,
      kappa = de_res$kappa,
      de_rd = de_rd_total
    )

    res_total <- list(
      IE_corrected = res_total_ie,
      DE_corrected = res_total_de
    )

    return(list(bias_only = res_bias_only, total = res_total))
  }


  # --- Run and Aggregate ---
  if (verbose) message(paste0("Starting ", B, " PBA iterations on ", n_cores, " core(s)..."))

  # # Select the correct iteration function based on the `bootstrap` flag
  # iteration_function <- if (bootstrap) {
  #   if (verbose) message("...using 2D Monte Carlo (bootstrap = TRUE).")
  #   pba_iter_bootstrap
  # } else {
  #   if (verbose) message("...using Normal Approximation (bootstrap = FALSE).")
  #   pba_iter_norm_approx
  # }

  # Run the single loop
  results_list <- if (n_cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(1:B, iteration_function, mc.cores = n_cores)
  } else {
    if (n_cores > 1 && .Platform$OS.type == "windows") {
      warning("mclapply is not supported on Windows. Using sequential 'lapply'.")
    }
    lapply(1:B, iteration_function)
  }

  # Separate the results
  bias_only_res <- lapply(results_list, `[[`, "bias_only")
  total_res <- lapply(results_list, `[[`, "total")

  # Filter out NULLs from failed iterations
  bias_only_res <- bias_only_res[!sapply(bias_only_res, is.null)]
  total_res <- total_res[!sapply(total_res, is.null)]

  if (length(bias_only_res) == 0 || length(total_res) == 0) {
    stop("All PBA iterations failed. Check inputs and prior functions.")
  }

  # Summarize both sets of results
  summary_bias_only <- summarize_pba_results(bias_only_res, "bias_only")
  summary_total <- summarize_pba_results(total_res, "total")

  # Combine into final data.tables
  final_ie_results <- rbind(summary_bias_only$IE_results, summary_total$IE_results)
  final_de_results <- rbind(summary_bias_only$DE_results, summary_total$DE_results)

  return(list(
    IE_results = final_ie_results,
    DE_results = final_de_results
  ))
}

