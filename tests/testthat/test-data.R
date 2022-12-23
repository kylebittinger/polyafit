test_that("Lung replicate samples data has documented types", {
  expect_equal(lapply(lungreplicate_samples, class), list(
    sample_id = "character", subject_id = "factor", replicate_id = "character",
    sample_type = "factor", study_group = "factor"))
})

test_that("Lung replicate samples data has documented dimensions", {
  expect_equal(dim(lungreplicate_samples), c(42, 5))
})

test_that("Lung replicate reads data has documented types", {
  expect_equal(lapply(lungreplicate_reads, class), list(
    otu_id = "factor", sample_id = "character", reads = "numeric",
    assignment = "character"))
})

test_that("Lung replicate reads data has documented dimensions", {
  expect_equal(dim(lungreplicate_reads), c(19362, 4))
})
