test_that("marginal probability of each category is beta binomial", {
  observed <- ppolya_marginal(c(1, 5, 6), c(8, 15, 16))
  expect_equal(observed[1], pbetabinom(1, 12, 8, 31)) # 15 + 16 = 31
  expect_equal(observed[2], pbetabinom(5, 12, 15, 24)) # 8 + 16 = 24
  expect_equal(observed[3], pbetabinom(6, 12, 16, 23)) # 8 + 15 = 23
})