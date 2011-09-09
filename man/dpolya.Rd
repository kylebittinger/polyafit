\name{dpolya}
\alias{dpolya}
\title{Compute multivariate Polya probabilities}
\usage{
  dpolya(x, alphas, log = TRUE)
}
\arguments{
  \item{x}{A vector of integers representing the number of
  times each outcome was observed.  The length must match
  that of alphas.}

  \item{alphas}{Parameters of the Dirichlet distribution}

  \item{log}{If TRUE, return the natural logarithm}
}
\value{
  Multivariate Polya probability mass function
}
\description{
  The multivariate Polya distribution, also known as the
  Dirichlet multinomial compound distribution, gives the
  probability of selecting \eqn{n_k} objects in each of
  \eqn{K} categories, when the probability of selecting an
  object from each category is given by a Dirichlet
  distribution with parameter \eqn{\alpha_k}: \deqn{
  P(\mathbf{n};\boldsymbol{\alpha}) =
  \frac{N!}{\prod_k{n_k!}}
  \frac{\Gamma(A)}{\prod_k{\Gamma(\alpha_k)}}
  \frac{\prod_k{\Gamma(\alpha_k + n_k)}}{\Gamma(A + N)}, }
  where \eqn{N = \sum_k n_k} and \eqn{A = \sum_k \alpha_k}.
}

