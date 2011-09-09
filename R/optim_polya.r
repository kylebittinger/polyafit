#' Estimate Dirichlet precision using Wicker's formula
#'
#' Returns an estimate of the Dirichlet precision (sum of the parameters of a
#' Dirichlet distribution) using a closed-form expression derived by Wicker
#' \emph{et. al.}
#'
#' @param props A matrix of training proportions, one column per category, one
#'   row per trial.
#' @param min_prop Minimum proportion used in the calculation.  Wicker's
#'   formula takes the log of each proportion, which fails if a proportion is
#'   zero.  This parameter provides a lower bound to proportions used in the
#'   formula.
#' @param min_value Minimum return value.  In some cases, Wicker's estimate
#'   is zero or negative.  This parameter is used to provide a lower bound to
#'   the return value.
#' @references
#'   Wicker, N., et al. "A Maximum Likelihood Approximation Method for 
#'   Dirichlet's Parameter Estimation," Comp. Stat. Data Anal. 52, 1315 (2008).
#' @return Estimated sum of parameters for Dirichlet distribution
#' @export
dirichlet_precision_wicker <- function(props, min_prop=1e-10, min_value=0.1) {
  if (is.vector(props)) {
    props <- matrix(props, nrow=1)
  }
  props[props < min_prop] <- min_prop
  n <- ncol(props)
  k <- nrow(props)
  f <- apply(props, 2, mean)
  euler_constant <- -digamma(1)
  
  term0 <- n * (k - 1) * euler_constant
  term1 <- n * sum(f * log(f))
  term2 <- sum(f * colSums(log(props)))
  wicker_estimate <- term0 / (term1 - term2)
  if (wicker_estimate < min_value) min_value else wicker_estimate
}


#' Find MLE of a multivariate Polya distribution
#'
#' @param counts A matrix of observed data, one column per category, one row
#'   per trial.
#' @return The result of the final call to optim().
#' @export
optim_polya <- function(counts) {

  polya_model <- function(alphas) {
    sum(apply(counts, 1, dpolya, alphas))
  }

  optim_step <- function(initial_params) {
    optim(
      initial_params,
      polya_model,
      method="L-BFGS-B",
      lower=rep(1e-10, length(initial_params)),
      upper=rep(1e10, length(initial_params)),
      control=list(maxit=1000, fnscale=-1, trace=0))
  }

  props <- t(apply(counts, 1, function (x) {x / sum(x)}))
  approx_precision <- dirichlet_precision_wicker(props, min_value=1)
  initial_params <- approx_precision * apply(props, 2, mean)

  results <- optim_step(initial_params)
  if (results$convergence != 0) {
    warning(
      "Convergence failed in call to optim(): ", 
      results$convergence,
      " ",
      results$message)
  }

  #
  # Check solution using variation of parameters
  #
  bad_variants_exist <- TRUE
  while (bad_variants_exist) {
    variant_multipliers <- expand.grid(
      idx  = seq_along(results$par), 
      mult = c(10, 2, 0.5, 0.1))
    variant_pars <- t(apply(variant_multipliers, 1, function (x) {
      idx <- x["idx"]
      mult <- x["mult"]
      results$par[idx] <- mult * results$par[idx]
    }))
    variant_likelihoods <- apply(variant_pars, 1, polya_model)
    bad_variants_exist <- max(variant_likelihoods) > results$value
    if (bad_variants_exist) {
      message("MLE exceeded with variation of parameters, re-optimizing...")
      top_par <- variant_pars[which.max(variant_likelihoods),]
      results <- optim_step(top_par)
    }
  }
  results
}
