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
