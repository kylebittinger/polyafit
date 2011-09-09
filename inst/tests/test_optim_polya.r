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
