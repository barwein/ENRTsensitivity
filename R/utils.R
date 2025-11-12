# Write note on this script -- utility function for the sensitivity analysis
# Should be a note documentation


### We focus only on identification of IE and DE and not the DE bounds!!! ###

# Libraries ---------------------------------------------------------------

# library(data.table)
# library(proxy)

# General functions -------------------------------------------------------

dist_lp_norm <- function(X_e, X_a = NULL, p){
  # L_p Norm distance

  # --- ONE MATRIX (output dim: n_e x n_e) ---
  if (is.null(X_a)) {
    method <- if (is.infinite(p)) "maximum" else
      if (p == 2) "euclidean" else
        if (p == 1) "manhattan" else "minkowski"
    return(as.matrix(stats::dist(X_e,
                                 method = method,
                                 p = if (method == "minkowski") p else 2)))
  }

  # --- TWO MATRICES (output dim: n_e x n_a) ---
  lp_fun <- function(a, b, p) {
    if (is.infinite(p)) max(abs(a - b)) else (sum(abs(a - b)^p))^(1 / p)
  }
  return(as.matrix(proxy::dist(X_e, X_a, method = lp_fun, p = p)))



  # if (is.null(X_a)) {X_a <- X_e}
  # EE <- rowSums(X_e^2) # length n_e
  # AA <- rowSums(X_a^2) # length n_a
  # DD_ <- -2 * X_e %*% t(X_a) + outer(EE, AA, `+`)
  # DD_[DD_ < 0] <- 0 # numerical stability
  # return(sqrt(DD_))
  # # return(as.matrix(dist(X)))
}

dist_inner_cosine <- function(X_e, X_a = NULL, dist){
  # Inner product: X_i^TX_j
  # And possible Cosine distance: 1 - (X_i^TX_j) / (||X_i|| ||X_j||)
  center <- FALSE
  if (is.null(X_a)) {
    X_a <- X_e
    center <- TRUE
  }
  inner_prod <- X_e %*% t(X_a)
  # Center the inner product matrix by removing the diagonal
  if (center){
    inner_prod_centered <- inner_prod - diag(diag(inner_prod))
  } else {
    inner_prod_centered <- inner_prod
  }
  # Return the appropriate distance metric
  if (dist == "inner"){
    return(inner_prod_centered)
  }
  if (dist == "cosine"){
    norms_e <- sqrt(rowSums(X_e^2)) # length n_e
    norms_a <- sqrt(rowSums(X_a^2)) # length n_a
    denom <- outer(norms_e, norms_a)
    # To avoid division by zero, set zero norms to a small value
    denom[denom == 0] <- .Machine$double.eps
    return(1 - inner_prod / denom)
    # return(1 - inner_prod_centered / denom)
    # return(inner_prod_centered / outer(norms, norms))
  }
  stop("Invalid distance metric specified.")
}

dist_func <- function(X_e, X_a = NULL, dist = "norm", p = 1){
  stopifnot(is.matrix(X_e))
  if (!is.null(X_a)) stopifnot(is.matrix(X_a), ncol(X_a) == ncol(X_e))
  if (!(is.infinite(p) || (is.numeric(p) && p >= 1)))
    stop("p must be >= 1 or Inf (unweighted).")

  if(dist == "norm"){
    return(dist_lp_norm(X_e, X_a, p))
  }
  if(dist %in% c("cosine","inner")){
    return(dist_inner_cosine(X_e, X_a, dist))
  }
  stop("Invalid distance metric specified.")
}


get_dist_matrix <- function(X_e,
                            X_a = NULL,
                            dist = "norm",
                            p = 1){

  if (!dist %in% c("norm","cosine","inner")){
    stop(paste0("Type=", dist, " of covariates distance is not supported."))
  }
  # Build the n_e x n_a distance matrix
  # Note that when X_a=NULL, this is the n_e x n_e distance matrix (with 0 diagonal)
  # and when X_a is given, this is the n_e x n_a distance matrix
  # with no restriction on the diagonal
  return(dist_func(X_e, X_a, dist, p))
}


# Weights and probs from distances ----------------------------------------------------

.col_logsumexp <- function(L){
  # Helper function to compute the column-wise LogSumExp of matrix L
  # L is either n_e x n_e or n_e x n_a
  # Were doing Softmax by column in L
  # For one column X: LSE(X) = log(sum(exp(X))) = max(X) + log(sum(exp(X - max(X))))
  # This is more numerically stable
  col_max <- apply(L, 2, max, na.rm=TRUE) # max by column
  Lc <- sweep(L, 2, col_max, "-") # X - max(X)
  Lc[!is.finite(Lc)] <- -Inf
  S  <- colSums(exp(Lc)) # sum(exp(X - max(X)))
  logsumexp <- col_max + log(S) # max(X) + log(sum(exp(X - max(X))))
  logsumexp[!is.finite(col_max)] <- -Inf
  return(logsumexp)
}

hetero_pi_weight_from_dist_ <- function(D,
                                        gamma = -1,
                                        ego_index = NULL
                                        ){
  # D: n_e x n_a (or n_e x n_e) distances matrix
  # gamma: scalar "temperature" (negative for distances; positive for inner prod)
  # If self_zero=TRUE and D is square, force P[ii] = 0 (excluded in the softmax)

  L <- gamma * D
  L[is.na(L)] <- -Inf # exclude missing distances
  # For ego-ego distance, make the diagonal -inf -> prob = 0
  if (!is.null(ego_index)){
    L[cbind(ego_index, seq_along(ego_index))] <- -Inf
  }
  if (is.null(ego_index) && nrow(D) == ncol(D)){
    diag(L) <- -Inf  # exclude self distances in ego-ego case
  }

  # Weights matrix: W_ij = exp(gamma * d_ij)
  # Required for Heterogeneous number of missing edges case
  weights <- exp(L)

  # Probs will be computed via Softmax by columns of W
  # We use LogSumExp trick for numerical stability
  # P_ij = exp(gamma*D_ij - LSE(D_j))
  lse <- .col_logsumexp(L)
  log_prob <- sweep(L, 2, lse, "-") # log-softmax by column
  log_prob[, is.infinite(lse)] <- -Inf # if all -Inf in col, set all probs to 0

  prob <- exp(log_prob)
  prob[!is.finite(prob)] <- 0

  return(list(weights = weights,
              prob = prob))
}

