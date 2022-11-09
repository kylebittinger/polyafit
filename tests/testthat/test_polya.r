test_that("equal to binomial coefficient when 2 parameters are given", {
  expect_equal(multichoose(c(1, 4)), lchoose(5, 4))
  expect_equal(multichoose(c(7, 3)), lchoose(10, 3))
  expect_equal(multichoose(c(3, 7)), lchoose(10, 3))
})

test_that("equal to N! when all parameters are 1", {
  expect_equal(multichoose(rep(1, 8)), lfactorial(8))
})


test_that("function is equal to lbeta when 2 parameters are given", {
  expect_equal(multibeta(c(5, 4)), lbeta(5, 4))
})

test_that("function is equal to 1/gamma(K) when all counts are 1", {
  expect_equal(multibeta(rep(1, 10)), -lgamma(10))
  expect_equal(multibeta(rep(1, 41)), -lgamma(41))
})

test_that("distribution is uniform when alphas are 1", {
  for (x in list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))) {
    expect_equal(dpolya(x, c(1, 1, 1)), log(1/3))
  }
})

test_that("distribution is beta-binomial when 2 parameters are given", {
  expect_equal(dpolya(c(1, 0), c(2, 1)), log(2/3))
  expect_equal(dpolya(c(0, 1), c(2, 1)), log(1/3))
})

test_that("distribution is sharply peaked when parameters are large", {
  expect_equal(dpolya(c(1,0,0), c(5000, 1, 1)), log(1), tolerance = 1e-3)
  expect_equal(dpolya(c(0,1,0), c(5000, 1, 1), log=FALSE), 0, tolerance = 1e-3)
})

test_that("distribution is normalized", {
  test_normalization <- function(perms, ntrials=10) {
    for (trial in 1:ntrials) {
      alphas <- runif(ncol(perms), min=0, max=5)
      probs <- apply(perms, 1, function (x) {
        dpolya(x, alphas, log=FALSE)
      })
      expect_equal(sum(probs), 1)
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
