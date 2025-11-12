
utils::globalVariables(c(
  ".", "ie_rd", "de_rd", "pi_param", "spec", "var_to_use",
  "ie_rd_boot_var", "ie_rd_var", "de_rd_boot_var", "de_rd_var",
  "std_err", "ci_low", "ci_high"
))

#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table :=
#' @importFrom data.table .N
#' @importFrom data.table .SD
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom proxy dist
#' @importFrom stats dist
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats qnorm
#' @importFrom stats rnorm
#' @importFrom stats predict
#' @importFrom stats quantile
#' @importFrom stats rpois
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom parallel mclapply
## usethis namespace: end
NULL
