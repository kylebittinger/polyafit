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

test_that("feature_enrichment gives correct result", {
  p <- pfit(x)
  pe <- feature_enrichment(p)
  expected_pe <- tibble::tibble(
    observation = c(
      "Sample A", "Sample A", "Sample A", "Sample A", "Sample A",
      "Sample B", "Sample B", "Sample B", "Sample B", "Sample B"),
    feature = c(
      "asv1", "asv2", "asv3", "asv4", "asv5",
      "asv1", "asv2", "asv3", "asv4", "asv5"),
    counts = c(10, 20, 60, 80, 120, 11, 25, 49, 300, 134),
    expected_counts = c(
      10.7515761195116, 19.6201486571511, 43.1461800195183,
      119.898315498794, 96.5837797050247, 19.2416138138845,
      35.113300527798, 77.2167842418277, 214.576640496118,
      172.851660920372),
    sigma = c(
      -0.086000458349907, 0.0326990083761397, 1.02392105818105,
      -1.75167260510861, 1.07418102470083, -0.543073663060615,
      -0.501340185350642, -0.987175411073878, 2.15970141483643,
      -1.02633457291013),
    p.value = c(
      0.412504752046599, 0.402301175300987, 0.147459692209273,
      0.961683706669587, 0.139743805230948, 0.628205565475373,
      0.631763754851917, 0.833639492127292, 0.0166259568998115,
      0.842497070103083))
  expect_equal(pe, expected_pe)
})
