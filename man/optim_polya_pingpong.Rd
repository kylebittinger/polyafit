\name{optim_polya_pingpong}
\alias{optim_polya_pingpong}
\title{Find MLE of a multivariate Polya distribution, optimizing mean and precision separately}
\usage{
  optim_polya_pingpong(counts, max_iter = 100)
}
\arguments{
  \item{counts}{A matrix of observed data, one column per
  category, one row per trial.}

  \item{max_iter}{Maximum number of iterations}
}
\value{
  The result of the final call to optim().
}
\description{
  Find MLE of a multivariate Polya distribution, optimizing
  mean and precision separately
}

