context("test-get_series")

extdata_dir <- system.file("extdata", package = "sigminer", mustWork = TRUE)

test_that("get_* series functions work", {
  skip_on_cran()

  cp <- read_copynumber(extdata_dir, pattern = "txt", genome_build = "hg19")
  expect_s4_class(cp, "CopyNumber")

  cn_list <- get_cnlist(cp)
  expect_true(inherits(cn_list, "list"))
  expect_true(length(cn_list) == 3)

  cn_features <- get_features(cn_list, cores = 1)
  expect_true(inherits(cn_features, "list"))
  expect_true(length(cn_features) == 6)

  cn_components <- get_components(cn_features)
  expect_true(inherits(cn_components, "list"))
  expect_true(length(cn_components) == 6)

  cn_matrix <- get_matrix(cn_features, cn_components)
  expect_equal(dim(cn_matrix) == c(3, 25), c(TRUE, TRUE))

  # parallel computation
  cn_features <- get_features(cn_list, cores = 2)
  expect_true(inherits(cn_features, "list"))
  expect_true(length(cn_features) == 6)
  cn_matrix2 <- get_matrix(cn_features, cn_components, cores = 2)
  expect_equal(dim(cn_matrix2) == c(3, 25), c(TRUE, TRUE))
})

extdata_dir <- NULL
