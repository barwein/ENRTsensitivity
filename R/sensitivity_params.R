# Source ------------------------------------------------------------------

# source("src/utils.R")

# Homogeneous case --------------------------------------------------------

#' @title Calculate Homogeneous Exposure Probabilities (pi)
#' @description
#' Calculates the probability of exposure to at least one treated ego (\eqn{\pi_i})
#' under a homogeneous contamination model. This model assumes that the
#' probability of a latent edge (or the number of latent edges) is uniform
#' across all possible ego-ego or alter-ego pairs.
#'
#' This function is used to generate sensitivity parameters for `enrt_sa` and
#' `enrt_pba` based on "Example 1 (Homogeneous probabilities)" and
#' "Example 2 (Homogeneous number of edges)" from the paper.
#'
#' @param rho_vec A numeric vector of homogeneous edge probabilities, \eqn{\rho^e} or
#'  \eqn{\rho^a}. Each value should be in [0, 1].
#' @param m_vec A numeric vector of the expected number of missing edges, \eqn{m^e} or
#'   \eqn{m^a}. Values must be non-negative.
#' @param n_e The total number of egos in the sample.
#' @param n_a The total number of alters in the sample (required if `type = "alter"`
#'   and `m_vec` is used).
#' @param type A string, either `"ego"` or `"alter"`.
#' * If `"ego"`, calculates \eqn{\pi_i^e = 1 - (1 - p_z \rho^e)^{n_e - 1}}.
#'  * If `"alter"`, calculates \eqn{\pi_i^a = p_z + (1 - p_z)[1 - (1 - p_z \rho^a)^{n_e - 1}]}.
#' @param pz The known probability of an ego being assigned to treatment, \eqn{Pr(Z=1)}.
#'
#' @return A named list where each element is a scalar \eqn{\pi} value.
#'   The names of the list correspond to the values from `rho_vec` or `m_vec`.
#'
#' @export
#'
#' @examples
#' # --- Example 1: Homogeneous Probabilities ---
#' # Assume 150 egos, pz = 0.5
#' # What is the exposure probability for egos if each ego-ego edge
#' # has a 0.5%, 1%, or 2% probability of existing?
#'
#' pi_e_list <- pi_homo(rho_vec = c(0.005, 0.01, 0.02),
#'                      n_e = 150,
#'                      type = "ego",
#'                      pz = 0.5)
#' print(pi_e_list)
#'
#' # --- Example 2: Homogeneous Number of Edges ---
#' # Assume 150 egos, 263 alters, pz = 0.5
#' # What is the exposure probability for alters if we expect
#' # m_a = 100, 200, or 300 total expected missing alter-ego edges?
#'
#' pi_a_list <- pi_homo(m_vec = c(100, 200, 300),
#'                      n_e = 150,
#'                      n_a = 263,
#'                      type = "alter",
#'                      pz = 0.5)
#' print(pi_a_list)
#'
pi_homo <- function(rho_vec = NULL,
                    m_vec = NULL,
                    n_e,
                    n_a = NULL,
                    type,
                    pz = 0.5
){

  # Assert that 'type' in c("ego","alter")
  if (!type %in% c("ego","alter")){
    stop("Invalid prob type; 'type' should be one of 'ego' or 'alter'")
  }
  # Edges prob with homogeneous prob or number of missing edges
  # If rho is a probability (Homogeneous probs case)
  # if(rho <= 1 && rho >= 0 && !is.integer(rho)){
  if(!is.null(rho_vec)){
    rho_ij_vec <- rho_vec
  } else if (!is.null(m_vec)){
    # If rho is an integer (Homogeneous number of missing edges case)
    if(any(m_vec < 0)){
      stop("All values in 'm_vec' must be non-negative.")
    }
    rho_ij_vec <- lapply(m_vec,
                         function(m){
                           ifelse(type == "ego",
                                  m / choose(n_e, 2),
                                  m / (n_a*(n_e-1)))
                         })
  } else{
    stop("Invalid inputs. One of 'rho_vec' or 'm_vec' should be supplied.")
  }
  # Prob unit is not exposed to additional treated ego
  # p_not_expos <- (1 - pz*rho_ij)^(n_e-1)
  # pi_ <- ifelse(type == "ego",
  #               1 - p_not_expos,
  #               pz + (1 - pz)*(1 - p_not_expos))
  # return(pi_)
  pi_list <- lapply(rho_ij_vec, function(rho_ij){
    p_not_expos <- (1 - pz*rho_ij)^(n_e-1)
    if (type == "ego"){
      1 - p_not_expos
    } else {
      pz + (1 - pz)*(1 - p_not_expos)
    }
  })
  if (!is.null(rho_vec)){
    # names(pi_list) <- paste0("rho=", rho_vec)
    names(pi_list) <- rho_vec
  } else {
    # names(pi_list) <- paste0("m=", m_vec)
    names(pi_list) <- m_vec
  }
  return(pi_list)
}


# Heterogeneous case ------------------------------------------------------

#' @title Calculate Heterogeneous Exposure Probabilities (pi)
#' @description
#' Calculates the probability of exposure to at least one treated ego (\eqn{\pi_i})
#' under a heterogeneous contamination model. This model assumes that the
#' probability of a latent edge depends on unit-level covariates \eqn{X},
#' typically through a distance or similarity function.
#'
#' This function generates sensitivity parameters for `enrt_sa` and
#' `enrt_pba` based on "Example 3 (Heterogeneous probabilities)" and
#' "Example 4 (Heterogeneous number of missing edges)" from the paper.
#'
#' @param X_e A numeric matrix of covariates for the egos (n_e rows).
#' @param X_a A numeric matrix of covariates for the alters (n_a rows).
#'   If `NULL`, the function calculates ego-ego probabilities (\eqn{\pi_i^e}).
#'   If provided, it calculates alter-ego probabilities (\eqn{\pi_i^a}).
#' @param m_vec A numeric vector of the expected number of missing edges, \eqn{m^e} or
#'   \eqn{m^a}. Used for "Example 4". If `NULL` (default),
#'   the function uses the "Example 3" model.
#' @param gamma A numeric vector of sensitivity parameters (\eqn{\gamma^e} or \eqn{\gamma^a})
#'   that control the influence of covariate similarity.
#'   If `m_vec` is provided, `gamma` must be a single scalar.
#' @param dist A string specifying the distance function.
#'   Supported: `"norm"` (Lp norm), `"cosine"`, `"inner"`.
#' @param ego_index A numeric vector mapping each alter to their ego's index.
#'   Required only when `type_ == "alter"` (i.e., `X_a` is not `NULL`).
#'   This is used to exclude distances to an alter's *own* ego.
#' @param pz The known probability of an ego being assigned to treatment, \eqn{Pr(Z=1)}.
#' @param p The power for the Lp norm (`dist = "norm"`), e.g., `p=2` for
#'   Euclidean (default) or `p=1` for Manhattan.
#' @param return_rho_ij Logical. If `TRUE` and `m_vec` is used, returns a
#'   list containing both the 'pi' vectors and the 'rho' (edge probability)
#'   matrices.
#'
#' @return
#' If `m_vec` is `NULL`: A named list where each element is a *vector*
#'   of \eqn{\pi_i} values (length n_e or n_a). The names correspond to `gamma` values.
#' If `m_vec` is not `NULL`: A named list where each element is a *vector*
#'   of \eqn{\pi_i} values. The names correspond to `m_vec` values.
#'
#' @importFrom stats dist
#' @importFrom proxy dist
#'
#' @export
#'
#' @examples
#' # --- Simulate data for examples ---
#' set.seed(123)
#' n_e <- 50
#' n_a <- 100
#' pz <- 0.5
#' X_e <- matrix(rnorm(n_e * 2), ncol = 2, dimnames = list(NULL, c("X1", "X2")))
#' X_a <- matrix(rnorm(n_a * 2), ncol = 2, dimnames = list(NULL, c("X1", "X2")))
#' ego_id_a <- sample(1:n_e, n_a, replace = TRUE) # Map alters to egos
#'
#' # --- Example 3: Heterogeneous Probabilities (Egos) ---
#' # Calculate pi_e vectors based on covariate similarity (gamma)
#'
#' pi_e_list_hetero <- pi_hetero(X_e = X_e,
#'                               gamma = c(0.5, 1.0, 2.0),
#'                               dist = "norm",
#'                               p = 2,
#'                               pz = pz)
#'
#' # Each element is a vector of length n_e
#' sapply(pi_e_list_hetero, length)
#' sapply(pi_e_list_hetero, mean)
#'
#' # --- Example 4: Heterogeneous Number of Edges (Alters) ---
#' # Calculate pi_a vectors based on an expected number of missing edges (m_a),
#' # where edge probability is weighted by covariate similarity (gamma = 1).
#'
#' pi_a_list_hetero <- pi_hetero(X_e = X_e,
#'                               X_a = X_a,
#'                               m_vec = c(50, 100),
#'                               gamma = 1.0,
#'                               dist = "norm",
#'                               ego_index = ego_id_a,
#'                               p = 2,
#'                               pz = pz)
#'
#' # Each element is a vector of length n_a
#' sapply(pi_a_list_hetero, length)
#' sapply(pi_a_list_hetero, mean)
#'
pi_hetero <- function(X_e,
                      X_a = NULL,
                      m_vec = NULL,
                      gamma = -1,
                      dist = "norm",
                      ego_index = NULL,
                      pz = 0.5,
                      p = 1,
                      return_rho_ij = FALSE){

  type_ <- ifelse(is.null(X_a), "ego", "alter")

  if (type_ == "alter" && ncol(X_e) != ncol(X_a)){
    stop("X_e and X_a must have the same number of columns.")
  }

  if (!is.null(m_vec) && length(gamma) > 1){
    stop("When 'm_vec' is given, 'gamma' must be a scalar.")
  }

  if (!is.null(m_vec) && any(m_vec < 0)){
    stop("All values in 'm_vec' must be non-negative.")
  }

  # Distance matrix
  D <- get_dist_matrix(X_e, X_a, dist, p)


  # get weights and probs matrices (for each gamma value)
  W_P_list <- lapply(gamma, function(g){
    hetero_pi_weight_from_dist_(D,
                                g,
                                ego_index)
  })

  if(is.null(m_vec)){
    # Heterogeneous probs case.
    pi_by_gamma <- lapply(W_P_list, function(wp){
      P <- wp$prob
      p_exposed <- 1 - apply(pz*P, 2, function(x) prod(1 - x))
      if (type_ == "ego"){
        p_exposed
      } else {
        pz + (1 - pz)*p_exposed
      }
    })
    # names(pi_by_gamma) <- paste0("gamma=", gamma)
    names(pi_by_gamma) <- gamma
    if (length(pi_by_gamma) == 1){
      return(pi_by_gamma[[1]])
    } else {
      return(pi_by_gamma)
    }
  } else{
    # Heterogeneous number of missing edges case
    # Compute all 'pi' probs by all 'm_vec' values at once
    W <- W_P_list[[1]]$weights
    W_sum <- ifelse(type_ == "ego",
                    sum(W[lower.tri(W)]),
                    sum(W))

    rho_by_m <- lapply(m_vec, function(m){m*W / W_sum})
    # names(rho_by_m) <- paste0("m=", m_vec)
    names(rho_by_m) <- m_vec

    # rho_ij <- missing_num*W / W_sum
    if (type_ == "ego"){
      # rho_ij <- rho_ij - diag(diag(rho_ij)) # force diagonal to 0
      rho_by_m <- lapply(rho_by_m, function(rho_mat){
        rho_mat - diag(diag(rho_mat))
      })
    }
    # Compute pi probs for each m
    pi_by_m <- lapply(rho_by_m, function(rho_ij){
      p_exposed <- 1 - apply(pz*rho_ij, 2, function(x) prod(1 - x))
      if (type_ == "ego"){
        p_exposed
      } else {
        pz + (1 - pz)*p_exposed
      }
    })
    if (return_rho_ij){
      return(list(pi = pi_by_m, rho = rho_by_m))
    } else{
      return(pi_by_m)
    }
  }
}


