\name{dirichlet_precision_wicker}
\alias{dirichlet_precision_wicker}
\title{Estimate Dirichlet precision using Wicker's formula}
\usage{
  dirichlet_precision_wicker(props, min_prop = 1e-10,
  min_value = 0.1)
}
\arguments{
  \item{props}{A matrix of training proportions, one column
  per category, one row per trial.}

  \item{min_prop}{Minimum proportion used in the
  calculation.  Wicker's formula takes the log of each
  proportion, which fails if a proportion is zero.  This
  parameter provides a lower bound to proportions used in
  the formula.}

  \item{min_value}{Minimum return value.  In some cases,
  Wicker's estimate is zero or negative.  This parameter is
  used to provide a lower bound to the return value.}
}
\value{
  Estimated sum of parameters for Dirichlet distribution
}
\description{
  Returns an estimate of the Dirichlet precision (sum of
  the parameters of a Dirichlet distribution) using a
  closed-form expression derived by Wicker \emph{et. al.}
}
\references{
  Wicker, N., et al. "A Maximum Likelihood Approximation
  Method for Dirichlet's Parameter Estimation," Comp. Stat.
  Data Anal. 52, 1315 (2008).
}

