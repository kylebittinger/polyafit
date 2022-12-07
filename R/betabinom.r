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

betabinom_pval <- function (k, n, a, b) {
  1 - pbetabinom(k, n, a, b, log.p = FALSE)
}

betabinom_mean <- function(n, a, b) {
  n * a / (a + b)
}

betabinom_variance <- function (n, a, b) {
  var_numerator <- n * a * b * (a + b + n)
  var_denominator <- ((a + b) ^ 2) * (a + b + 1)
  var_numerator / var_denominator
}

betabinom_sd <- function (n, a, b) {
  sqrt(betabinom_variance(n, a, b))
}

betabinom_expected <- function (k, n, a, b) {
  betabinom_mean(n, a, b)
}

betabinom_sigma <- function (k, n, a, b) {
  bbm <- betabinom_mean(n, a, b)
  bbsd <- betabinom_sd(n, a, b)
  (k - bbm) / bbsd
}
