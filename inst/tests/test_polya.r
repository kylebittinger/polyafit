context("multichoose")

test_that("equal to binomial coefficient when 2 parameters are given", {
  expect_that(multichoose(c(1, 4)), equals(lchoose(5, 4)))
  expect_that(multichoose(c(7, 3)), equals(lchoose(10, 3)))
  expect_that(multichoose(c(3, 7)), equals(lchoose(10, 3)))
})

test_that("equal to N! when all parameters are 1", {
  expect_that(multichoose(rep(1, 8)), equals(lfactorial(8)))
})


context("multibeta")

test_that("function is equal to lbeta when 2 parameters are given", {
  expect_that(multibeta(c(5, 4)), equals(lbeta(5, 4)))
})

test_that("function is equal to 1/gamma(K) when all counts are 1", {
  expect_that(multibeta(rep(1, 10)), equals(-lgamma(10)))
  expect_that(multibeta(rep(1, 41)), equals(-lgamma(41)))
})


context("dpolya")

test_that("distribution is uniform when alphas are 1", {
  for (x in list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))) {
    expect_that(dpolya(x, c(1, 1, 1)), equals(log(1/3)))
  }
})

test_that("distribution is beta-binomial when 2 parameters are given", {
  expect_that(dpolya(c(1, 0), c(2, 1)), equals(log(2/3)))
  expect_that(dpolya(c(0, 1), c(2, 1)), equals(log(1/3)))
})

test_that("distribution is sharply peaked when parameters are large", {
  expect_that(dpolya(c(1,0,0), c(5000, 1, 1)), equals(log(1), tol=1e-3))
  expect_that(dpolya(c(0,1,0), c(5000, 1, 1), log=FALSE), equals(0, tol=1e-3))
})

test_that("distribution is normalized", {
  test_normalization <- function(perms, ntrials=10) {
    for (trial in 1:ntrials) {
      alphas <- runif(ncol(perms), min=0, max=5)
      probs <- apply(perms, 1, function (x) {
        dpolya(x, alphas, log=FALSE)
      })
      expect_that(sum(probs), equals(1))
    }
  }
  # N = 2, K = 3
  test_normalization(matrix(
    c(2, 0, 0, 
      1, 1, 0, 
      1, 0, 1, 
      0, 2, 0, 
      0, 1, 1, 
      0, 0, 2),
    ncol=3, byrow=TRUE))
  # N = 1, K = 3
  test_normalization(matrix(
    c(1, 0, 0, 
      0, 1, 0, 
      0, 0, 1),
    ncol=3, byrow=TRUE))
  # N = 5, K = 2
  test_normalization(matrix(
    c(0, 5,
      1, 4,
      2, 3,
      3, 2,
      4, 1,
      5, 0),
    ncol=2, byrow=TRUE))
  # N = 4, K = 3
  test_normalization(matrix(
    c(4, 0, 0, 
      3, 1, 0, 
      3, 0, 1, 
      2, 2, 0, 
      2, 1, 1, 
      2, 0, 2, 
      1, 3, 0, 
      1, 2, 1, 
      1, 1, 2, 
      1, 0, 3, 
      0, 4, 0, 
      0, 3, 1, 
      0, 2, 2, 
      0, 1, 3, 
      0, 0, 4),
      ncol=3, byrow=TRUE))
})
