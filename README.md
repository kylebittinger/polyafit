
<!-- README.md is generated from README.Rmd. Please edit that file -->

# polyafit

<!-- badges: start -->
<!-- badges: end -->

`polyafit` is an R package that provides tools to fit count data to a
[Dirichlet-multinomial
distribution](https://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution).
This distribution is sometimes called the multivariate Pólya
distribution, from which the package name is derived. The package was
created specifically for 16S rRNA marker gene sequence data from the
human microbiome, but in principle could be applied to other types of
data.

## Installation

You can install polyafit from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kylebittinger/polyafit")
```

## Example

``` r
library(tidyverse)
```

We will simulate two samples, and check to see if any features are
different, based on the distribution of all the features in the samples.
We assume that most of the features come from the same underlying
distribution, so we can use the full set of features to estimate the
amount that feature proportions vary between samples.

First, we’ll generate data for both samples using the same underlying
distribution. We will randomly generate some proportions using the
Chi-Squared distribution. The numbers won’t add to 1, so we’ll use a
convenience function to divide by the sum.

``` r
as_proportion <- function (x) x / sum(x)
set.seed(1)
example_props <- as_proportion(rchisq(50, 2))
```

Now, we’ll use a multinomial distribution to draw two samples with the
same underlying proportions. We’ll call these samples lung and OW, to
represent lung (lung) and oral (oral) microbiome samples from the same
individual.

Just to demonstrate that our method doesn’t require the same number of
counts in each sample, we’ll draw a different number of items for each
sample.

``` r
lung_cts <- rmultinom(1, 500, example_props)
oral_cts <- rmultinom(1, 600, example_props)
```

Finally, we’ll remove features with fewer than 5 counts total. This
method does not work well for features that appear sparsely in the
sample set.

``` r
keep <- (lung_cts + oral_cts) >= 5
lung_cts <- lung_cts[keep]
oral_cts <- oral_cts[keep]
```

### Analysis: no outlier

Run an outlier analysis on the counts. Here, the features all come from
the same underlying distribution, and we expect to see no outliers.

``` r
tibble(Lung = lung_cts, Oral = oral_cts) %>%
  mutate(across(everything(), as_proportion)) %>%
  ggplot(aes(y=Lung, x=Oral)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal()
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

Raw p-values for each feature. All p-values are &gt; 0.05.

``` r
pfit <- optim_polya(rbind(lung_cts, oral_cts))
pvals <- 1 - ppolya_marginal(lung_cts, pfit$par, log.p = FALSE)
pvals
#>  [1] 0.57594911 0.29223480 0.25754418 0.75969231 0.75420050 0.45033511
#>  [7] 0.70868359 0.12255148 0.47080502 0.50268617 0.36401274 0.64525371
#> [13] 0.30713523 0.19694488 0.50823963 0.61460214 0.19874754 0.14448560
#> [19] 0.39422966 0.39828908 0.25869333 0.08257417 0.33150315 0.25598698
#> [25] 0.49767320 0.66720092 0.30311445 0.26502399 0.25754418 0.09918343
#> [31] 0.44970098 0.39828908 0.40133580 0.31393493 0.88334831 0.93138181
#> [37] 0.27694878 0.49306991 0.47080502
```

The fitted object, `pfit`, has a vector of estimated parameter values
for the Dirichlet-multinomial (DM) distribution, which serves as the
mathematical model for our data. From the estimated parameter values, we
can compute an overdispersion parameter, theta. Theta tells us how much
“extra variance” the estimated distribution has, beyond what we expect
from a multinomial. Because the data were generated from a perfect
multinomial distribution, the estimated overdispersion parameter should
be very small.

``` r
1 / sum(pfit$par)
#> [1] 1.387863e-06
```

### Analysis: one outlier

Adjust one of the counts to be different in lung. Put in 150 counts for
the 4th feature.

``` r
lung_cts[4] <- 150
```

``` r
tibble(Lung = lung_cts, Oral = oral_cts) %>%
  mutate(across(everything(), as_proportion)) %>%
  ggplot(aes(y=Lung, x=Oral)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal()
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

Doing the same analysis, we can see that the p-value for the 4th feature
is very very low (&lt;&lt; 0.05).

``` r
pfit <- optim_polya(rbind(lung_cts, oral_cts))
pvals <- 1 - ppolya_marginal(lung_cts, pfit$par, log.p = FALSE)
pvals
#>  [1] 7.222374e-01 5.207786e-01 4.444319e-01 5.757440e-09 8.076512e-01
#>  [6] 5.795736e-01 7.079344e-01 2.759365e-01 5.808243e-01 6.724547e-01
#> [11] 5.387818e-01 7.240702e-01 5.025280e-01 3.742216e-01 6.508331e-01
#> [16] 7.455275e-01 3.499326e-01 2.932627e-01 5.204756e-01 5.157665e-01
#> [21] 4.506156e-01 2.017031e-01 5.532490e-01 4.843211e-01 5.838542e-01
#> [26] 6.622745e-01 5.084936e-01 6.106616e-01 4.444319e-01 3.362911e-01
#> [31] 6.662710e-01 5.157665e-01 5.087697e-01 4.767358e-01 8.426012e-01
#> [36] 9.376779e-01 5.757478e-01 6.211877e-01 5.808243e-01
```

Now, we re-calculate the overdispersion parameter, theta. It has
increased by three orders of magnitude.

``` r
1 / sum(pfit$par)
#> [1] 0.001526456
```
