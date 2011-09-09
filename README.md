# polyafit

`polyafit` is an R package that provides tools to fit data to a [multivariate
Polya distribution](http://en.wikipedia.org/wiki/Multivariate_Polya_distribution).

## Example

The input data is expected to be a matrix one row per trial, and one column per
category.

    > a <- matrix(c(5, 4, 7, 10, 12, 15, 4, 6, 2, 39, 43, 51), nrow=2)

The maximum likelihood estimate of the multivariate Polya parameters is returned
by polya_optim.

    > afit <- polya_optim(a)
    > afit$par
    [1]  2.363477  4.033356  6.164364  2.555273  5.856060 20.286586
