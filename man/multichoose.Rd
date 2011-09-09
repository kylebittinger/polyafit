\name{multichoose}
\alias{multichoose}
\title{Multinomial coefficient}
\usage{
  multichoose(x, log = TRUE)
}
\arguments{
  \item{x}{Vector representing the number of times each
  outcome was observed}

  \item{log}{If TRUE, return the natural logarithm}
}
\value{
  The multinomial coefficient
}
\description{
  Multinomial coefficients give the number of ways to
  partition a set of \eqn{N} objects into \eqn{K}
  containers, where each contains \eqn{n_k} elements.
  \deqn{{N \choose n_1, n_2, \ldots, n_K} =
  \frac{N!}{n_1!n_2! \cdots n_K!}}
}

