# rdirichlet copied from gregtools, Original code posted by Ben Bolker to
# R-News on Fri Dec 15 2000 [1].  Ben attributed the code to Ian Wilson
# i.wilson@maths.abdn.ac.uk. Subsequent modifications by Gregory R. Warnes
# greg@warnes.net.
# [1] http://www.r-project.org/nocvs/mail/r-help/2000/3865.html. 
rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol=l, byrow=TRUE)
  sm <- x %*% rep(1, l)
  x / as.vector(sm)
}


context("dirichlet_precision_wicker")

test_that("estimated precision is close to results in paper, Fig. 1", {
  alphas <- c(0.05, 0.1, 0.15, 0.2)
  true_precision <- sum(alphas)
  props <- rdirichlet(100, alphas)
  est_precision <- dirichlet_precision_wicker(props)
  # median and tolerance estimated from paper
  expect_that(est_precision - true_precision, equals(-0.3, tol=0.5))
})

test_that("dirichlet precision estimate works for single vector" ,{
  props <- c(0.1, 0.2, 0.6, 0.1)
  est_precision <- dirichlet_precision_wicker(props, min_value=-1e10)
  expect_that(est_precision, equals(0, tol=0.01))
})


context("optim_polya")

test_that("optimized parameters are good fit to those of input data", {
  alphas <- c(1.2, 3.3, 10.1, 0.5, 2.7)
  probs <- rdirichlet(50, alphas)
  counts <- t(apply(probs, 1, function (x) {rmultinom(1, 300, x)}))
  optim_result <- optim_polya(counts)
  # Don't constrain intercept -- remain more sensitive to incorrect results
  optim_fit <- summary(lm(optim_result$par ~ alphas))
  optim_slope <- optim_fit$coefficients[2, 1]
  optim_error <- optim_fit$coefficients[2, 2]
  optim_pval <- optim_fit$coefficients[2, 4]
  optim_intercept <- optim_fit$coefficients[1, 1]
  expect_that(optim_slope, equals(1, tol=0.5))
  expect_that(optim_error, equals(0, tol=0.1))
  expect_that(optim_pval, equals(0, tol=0.05))
  expect_that(optim_intercept, equals(0, tol=0.5))
})

test_that("true max. likelihood is found for problematic dataset", {
  counts <- matrix(c(
    73, 9, 23, 25, 82, 6, 2, 6, 534, 5, 1211, 13, 27, 8, 3, 950, 29, 11, 18, 
    50, 23, 6, 54, 274, 66, 27, 15, 24, 10, 13, 59, 0, 8, 630, 10, 6, 35, 22, 
    18, 740, 5, 54, 10, 11, 59, 14, 33, 1224, 28, 304, 11, 9, 14, 14, 4, 309, 
    5, 2058, 418, 21, 15, 18, 25, 54, 5, 108, 8, 80, 6, 304, 330, 7, 10, 23, 
    126, 8, 3, 1700, 250, 4, 7, 493, 7, 66, 122, 17, 729, 306, 157, 49, 68, 7, 
    448, 57, 917, 5, 37, 96, 32, 175, 111, 13, 3, 11, 128, 5, 35, 18, 68, 29, 
    310, 7, 3, 6, 3, 34, 0, 0, 4, 4, 231, 7, 578, 23, 2, 10, 3, 37, 28, 0, 16, 
    2, 21, 0, 14, 3, 29, 1, 10, 0, 5, 8, 0, 6, 0, 217, 1, 2, 21, 6, 65, 269, 1, 
    39, 11, 0, 26, 1, 10, 33, 2, 87, 16, 1, 0, 0, 4, 14, 1, 315, 0, 1, 0, 12, 
    35, 29, 5, 7, 1, 0, 1, 225, 156, 0, 16, 1, 1, 11, 8, 372, 629, 7, 1, 0, 1, 
    3, 24, 0, 85, 0, 50, 6, 56, 1, 3, 5, 39, 1, 7, 97, 0, 1, 28, 1, 3, 9, 4, 2, 
    3, 1, 34, 6, 72, 4), nrow=2, byrow=TRUE)
  op <- optim_polya(counts)
  # single call to optim() will find max likelihood of -1024.34
  expect_that(op$value, equals(-936.1, tol=1))
})

test_that("proportions are in reasonable order for problematic dataset", {
  counts <- matrix(c(
    45, 23, 13, 3, 7, 9, 7, 0, 13, 6, 8, 2, 110, 9, 32, 
    8, 92, 64, 73, 2, 31, 2, 6, 1, 19, 6, 196, 14, 415, 63, 61, 7, 
    2, 5, 4, 3, 661, 214, 22, 2, 124, 61, 4, 10, 37, 16, 438, 457, 
    59, 4, 1383, 935, 83, 118, 35, 25, 5, 2, 3, 6, 142, 9, 51, 8, 
    737, 723, 92, 12, 47, 19, 206, 126, 444, 149, 2, 7, 54, 4, 41, 
    11, 14, 9, 4, 2, 69, 5, 4, 3, 3, 10, 27, 8, 8, 0, 6, 0, 3, 3, 
    130, 11, 8, 2, 8, 1, 524, 134, 19, 26, 1, 6, 7, 1, 6, 7, 64, 
    41, 247, 179, 49, 28, 22, 1, 7, 4, 16, 7, 221, 90, 55, 12, 32, 
    27, 8, 3, 68, 1, 13, 1, 216, 181, 53, 12, 16, 1, 6, 0, 50, 30, 
    2, 7, 5, 2, 18, 14, 1260, 923, 11, 6, 4, 2, 8, 1, 1, 6, 72, 54, 
    7, 0, 107, 5, 34, 10, 89, 106, 245, 98, 4, 2, 98, 52, 17, 40, 
    22, 3, 50, 4, 131, 6, 78, 6, 1292, 714, 150, 5, 86, 6, 1862, 
    230, 594, 53, 154, 76, 693, 473, 1988, 224, 64, 0, 11, 7, 185, 
    50, 46, 8, 66, 89, 140, 103, 64, 52, 8, 1, 10, 3, 9, 2, 25, 40, 
    15, 5, 19, 3, 720, 469, 247, 98, 8, 0, 641, 70, 1439, 1074, 40, 
    1, 10, 1, 50, 2, 29, 0), nrow=2)
  par <- optim_polya(counts)$par

  # Columns with zeros in second row
  # Col   4 (x =  7)
  # Col  47 (x =  8)
  # Col  48 (x =  6)
  # Col  73 (x =  6)
  # Col  84 (x =  7)
  # Col 104 (x = 64)
  # Col 119 (x =  8)
  # Col 125 (x = 29)
  
  expect_that(par[104] > par[125], is_true()) # x=64 > x=29 
  expect_that(par[125] > par[ 47], is_true()) # x=29 > x=8
  expect_that(par[ 47] > par[  4], is_true()) # x=8 > x=7
  expect_that(par[  4] > par[ 48], is_true()) # x=7 > x=6
  expect_that(par[ 48], equals(par[ 73])) # x=6 matches
  expect_that(par[  4], equals(par[ 84])) # x=7 matches
  expect_that(par[ 47], equals(par[119])) # x=8 matches
})


test_that("repeated columns of counts have equal parameter estimates", {
  counts <- matrix(
    c(5, 20, 78, 1, 20,
      8,  4, 65, 3,  4),
    nrow=2, byrow=TRUE)
  op <- optim_polya(counts)
  # columns 2 and 5 are identical
  expect_that(op$par[2], equals(op$par[5], tol=0.01))
})

test_that("optimization works for single row of observations", {
  counts <- matrix(rep(1, 6), nrow=1)
  op <- optim_polya(counts)
  # Large tolerance essentially tests for upper bound of 19,999
  expect_that(op$par[1], equals(10000, tol=9999))
})
