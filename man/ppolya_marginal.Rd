% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polya_marginal.r
\name{ppolya_marginal}
\alias{ppolya_marginal}
\title{Compute the marginal CDF for each component of a multivariate Polya distribution}
\usage{
ppolya_marginal(x, alphas, log.p = TRUE)
}
\arguments{
\item{x}{Vector of observed counts in each category.}

\item{alphas}{Parameters of the distribution.}

\item{log.p}{If TRUE, return the natural logarithm.}
}
\value{
The marginal cumulative probability of observing the given number of
  counts for each component.
}
\description{
Compute the marginal CDF for each component of a multivariate Polya distribution
}
