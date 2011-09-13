\name{dbetabinom}
\alias{dbetabinom}
\title{Beta-binomial distribution: probability mass function}
\usage{
  dbetabinom(k, n, a, b, log = TRUE)
}
\arguments{
  \item{k}{Vector of successes.}

  \item{n}{Number of trials.}

  \item{a}{First shape parameter of Beta distribution (must
  be positive)}

  \item{b}{Second shape parameter of Beta distribution
  (must be positive)}

  \item{log}{If TRUE, return the natural logarithm}
}
\value{
  The probability of observing k successes in n trials
}
\description{
  Beta-binomial distribution: probability mass function
}

