# polyafit

`polyafit` is an R package that provides tools to fit data to a [multivariate
Polya distribution](http://en.wikipedia.org/wiki/Multivariate_Polya_distribution).

## Example

The input data is expected to be a matrix one row per trial, and one column per
category.

    > a <- matrix(c(5, 2, 8, 10, 19, 15, 4, 7, 0, 39, 43, 81), nrow=3)

The maximum likelihood estimate of the multivariate Polya parameters is returned
by polya_optim.

    > b <- polya_optim(a)
    > b$par
    [1]  3.878663 11.057931  2.507630 39.071530
