test_that("pfit makes an object from a matrix", {
  x <- matrix(c(1,2,3,4,5,6), nrow = 2)
  p <- pfit(x)
  expect_equal(p$data, x)
  expect_equal(class(p), "pfit")
})

test_that("pfit makes an object from a data frame", {
  x <- matrix(c(1,2,3,4,5,6), nrow = 2)
  rownames(x) <- c("A", "B")
  colnames(x) <- c("f1", "f2", "f3")
  xdf <- data.frame(
    sample_id = c("A", "B"), f1 = c(1, 2), f2 = c(3, 4), f3 = c(5, 6))
  p <- pfit(xdf)
  expect_equal(p$data, x)
  expect_equal(class(p), "pfit")
})

test_that("pfit_enrichment gives p-values", {
  x <- matrix(
    c(10, 11, 20, 25, 60, 49, 80, 300, 100, 114), nrow = 2,
    dimnames = list(c("Sample A", "Sample B"), c("f1", "f2", "f3", "f4", "f5")))
  p <- pfit(x)
  pe <- pfit_enrichment(p)
  expected_pe <- tibble::tibble(
    observation_idx = rep(1:2, each = 5),
    feature_idx = rep(1:5, times = 2),
    p_value = c(
    0.400446058361298, 0.387080796289815, 0.128409562332575, 0.961629256340697,
    0.15615678254315, 0.635050326369183, 0.642147372530849, 0.847828699271596,
    0.0166951712683215, 0.82804723768026))
  expect_equal(pe, expected_pe)
})
