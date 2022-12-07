#' Compute the marginal CDF for each component of a multivariate Polya distribution
#'
#' @param x Vector of observed counts in each category.
#' @param alphas Parameters of the distribution.
#' @param log.p If TRUE, return the natural logarithm.
#' @return The marginal cumulative probability of observing the given number of
#'   counts for each component.
#' @export
ppolya_marginal <- function(x, alphas, log.p=TRUE) {
  ppolya_apply(x, alphas, pbetabinom, log.p=log.p)
}

# Apply a function to each component of a Dirichlet-multinomial distribution.
# The function should have the signature f(k, n, a, b, ...), where k is the
# observed counts for the component.
ppolya_apply <- function(x, alphas, fcn, ...) {
  sum_a <- sum(alphas)
  sum_x <- sum(x)
  sapply(seq_along(x), function (idx) {
    x1 <- x[idx]
    a1 <- alphas[idx]
    fcn(x1, sum_x, a1, sum_a - a1, ...)
  })
}
