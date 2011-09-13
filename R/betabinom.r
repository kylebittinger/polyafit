#' Beta-binomial distribution: probability mass function
#'
#' @param k Vector of successes.
#' @param n Number of trials.
#' @param a First shape parameter of Beta distribution (must be positive)
#' @param b Second shape parameter of Beta distribution (must be positive)
#' @param log If TRUE, return the natural logarithm
#' @return The probability of observing k successes in n trials
#' @export
dbetabinom <- function(k, n, a, b, log=TRUE) {
  ldbb <- lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b)
  if (log) ldbb else exp(ldbb)
}

#' Beta-binomial distribution: cumulative distribution function
#'
#' @param k Vector of successes.
#' @param n Number of trials.
#' @param a First shape parameter of Beta distribution (must be positive)
#' @param b Second shape parameter of Beta distribution (must be positive)
#' @param log.p If TRUE, return the natural logarithm
#' @return The cumulative probability of observing k successes in n trials
#' @export
pbetabinom <- function(k, n, a, b, log.p=TRUE) {
  pbb <- sum(sapply(k, function (x) dbetabinom(0:x, n, a, b, log=FALSE)))
  if (log.p) log(pbb) else pbb
}