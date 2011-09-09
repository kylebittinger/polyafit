#' Multinomial coefficient
#'
#' Multinomial coefficients give the number of ways to partition a set of
#' \eqn{N} objects into \eqn{K} containers, where each contains \eqn{n_k} 
#' elements.
#' \deqn{{N \choose n_1, n_2, \ldots, n_K} = \frac{N!}{n_1!n_2! \cdots n_K!}}
#'
#' @param x Vector representing the number of times each outcome was observed
#' @param log If TRUE, return the natural logarithm
#' @return The multinomial coefficient
#' @export
multichoose <- function(x, log = TRUE) {
  lmc <- lfactorial(sum(x)) - sum(lfactorial(x)) 
  if (log) lmc else exp(lmc)
}

#' Multinomial Beta function
#' 
#' The multinomial extension of the Beta function is
#' \deqn{
#'   \mathrm{B}(\boldsymbol{\alpha}) = 
#'   \frac{\prod_i \Gamma(\alpha_i)}{\Gamma(\sum_i \alpha_i)}
#' }
#'
#' @param alphas Parameters of the multivariate Beta function
#' @param log If TRUE, return the natural logarithm
#' @return The multinomial Beta function
#' @export
multibeta <- function(alphas, log = TRUE) {
  lmb <- sum(lgamma(alphas)) - lgamma(sum(alphas))
  if (log) lmb else exp(lmb)
}

#' Compute multivariate Polya probabilities
#'
#' The multivariate Polya distribution, also known as the Dirichlet
#' multinomial compound distribution, gives the probability of selecting
#' \eqn{n_k} objects in each of \eqn{K} categories, when the probability of 
#' selecting an object from each category is given by a Dirichlet distribution
#' with parameter \eqn{\alpha_k}:
#' \deqn{
#'   P(\mathbf{n};\boldsymbol{\alpha}) = 
#'   \frac{N!}{\prod_k{n_k!}}
#'   \frac{\Gamma(A)}{\prod_k{\Gamma(\alpha_k)}}
#'   \frac{\prod_k{\Gamma(\alpha_k + n_k)}}{\Gamma(A + N)},
#' }
#' where \eqn{N = \sum_k n_k} and \eqn{A = \sum_k \alpha_k}.
#'
#' @param x A vector of integers representing the number of times each 
#'   outcome was observed.  The length must match that of alphas.
#' @param alphas Parameters of the Dirichlet distribution
#' @param log If TRUE, return the natural logarithm
#' @return Multivariate Polya probability mass function
#' @export
dpolya <- function(x, alphas, log = TRUE) {
  ldp <- multichoose(x) + multibeta(alphas + x) - multibeta(alphas)
  if (log) ldp else exp(ldp)
}
