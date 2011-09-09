\name{multibeta}
\alias{multibeta}
\title{Multinomial Beta function}
\usage{
  multibeta(alphas, log = TRUE)
}
\arguments{
  \item{alphas}{Parameters of the multivariate Beta
  function}

  \item{log}{If TRUE, return the natural logarithm}
}
\value{
  The multinomial Beta function
}
\description{
  The multinomial extension of the Beta function is \deqn{
  \mathrm{B}(\boldsymbol{\alpha}) = \frac{\prod_i
  \Gamma(\alpha_i)}{\Gamma(\sum_i \alpha_i)} }
}

