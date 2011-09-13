#' Compute the marginal CDF for each component of a multivariate Polya distribution
#'
#' @param x Vector of observed counts in each category.
#' @param alphas Parameters of the distribution.
#' @param log.p If TRUE, return the natural logarithm.
#' @return The marginal cumulative probability of observing the given number of
#'   counts for each component.
#' @export
ppolya_marginal <- function(x, alphas, log.p=TRUE) {
  sum_a <- sum(alphas)
  sum_x <- sum(x)
  sapply(seq_along(x), function (idx) {
    x1 <- x[idx]
    a1 <- alphas[idx]
    pbetabinom(x1, sum_x, a1, sum_a - a1, log.p=log.p)
  })
}
