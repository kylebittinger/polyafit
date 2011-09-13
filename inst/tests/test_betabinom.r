context("dbetabinom")

test_that("is normalized", {
  expect_that(sum(dbetabinom(0:10, 10, 4, 16, log=FALSE)), equals(1))
})

test_that("is zero when k > n", {
  expect_that(dbetabinom(7, 5, 3, 8, log=FALSE), equals(0))
})


context("pbetabinom")

test_that("is unity when k = n", {
  expect_that(pbetabinom(6, 6, 2, 34, log=FALSE), equals(1))
})

test_that("is sum of discrete probabilities", {
  expect_that(
    pbetabinom(2, 6, 9, 4, log=FALSE), 
    equals(sum(dbetabinom(c(0, 1, 2), 6, 9, 4, log=FALSE))))
})

test_that("is less than unity by sum of discrete probabilities", {
  expect_that(pbetabinom(4, 6, 9, 4, log=FALSE), 
    equals(1 - sum(dbetabinom(c(5, 6), 6, 9, 4, log=FALSE))))  
})
