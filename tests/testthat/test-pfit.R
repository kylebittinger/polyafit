x <- matrix(
  c(10, 11, 20, 25, 60, 49, 80, 300, 120, 134),
  nrow = 2,
  dimnames = list(
    c("Sample A", "Sample B"),
    c("asv1", "asv2", "asv3", "asv4", "asv5")))

xdf <- tibble::tibble(
  sample_id = c("Sample A", "Sample B"),
  asv1 = c(10, 11),
  asv2 = c(20, 25),
  asv3 = c(60, 49),
  asv4 = c(80, 300),
  asv5 = c(120, 134))

test_that("pfit makes an object from a matrix", {
  p <- pfit(x)
  expect_equal(p$data, x)
  expect_equal(class(p), "pfit")
})

test_that("pfit makes an object from a data frame", {
  p <- pfit(xdf)
  expect_equal(p$data, x)
  expect_equal(class(p), "pfit")
})

test_that("feature_enrichment gives p-values", {
  p <- pfit(x)
  pe <- feature_enrichment(p)
  expected_pe <- tibble::tibble(
    observation = c(
      "Sample A", "Sample A", "Sample A", "Sample A", "Sample A",
      "Sample B", "Sample B", "Sample B", "Sample B", "Sample B"),
    feature = c(
      "asv1", "asv2", "asv3", "asv4", "asv5",
      "asv1", "asv2", "asv3", "asv4", "asv5"),
    p.value = c(
      0.412504752046599, 0.402301175300987, 0.147459692209273,
      0.961683706669587, 0.139743805230948, 0.628205565475373,
      0.631763754851917, 0.833639492127292, 0.0166259568998115,
      0.842497070103083))
  expect_equal(pe, expected_pe)
})
