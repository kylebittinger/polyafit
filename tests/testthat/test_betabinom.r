test_that("is normalized", {
  expect_equal(sum(dbetabinom(0:10, 10, 4, 16, log=FALSE)), 1)
})

test_that("is zero when k > n", {
  expect_equal(dbetabinom(7, 5, 3, 8, log=FALSE), 0)
})


test_that("is unity when k = n", {
  expect_equal(pbetabinom(6, 6, 2, 34, log=FALSE), 1)
})

test_that("is sum of discrete probabilities", {
  expect_equal(
    pbetabinom(2, 6, 9, 4, log=FALSE),
    sum(dbetabinom(c(0, 1, 2), 6, 9, 4, log=FALSE)))
})

test_that("is less than unity by sum of discrete probabilities", {
  expect_equal(pbetabinom(4, 6, 9, 4, log=FALSE),
    1 - sum(dbetabinom(c(5, 6), 6, 9, 4, log=FALSE)))
})
