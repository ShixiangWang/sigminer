# Test for issue: NULL row names in NMF matrix when using sig_tally() with Wang method
# https://github.com/ShixiangWang/sigminer/issues/XXX

test_that("sig_tally with Wang method preserves row names for single sample", {
  # Load test data
  load(system.file("extdata", "toy_segTab.RData",
                   package = "sigminer", mustWork = TRUE))
  
  # Set seed for reproducibility
  set.seed(1234)
  
  # Add minor_cn column
  segTabs$minor_cn <- sample(c(0, 1), size = nrow(segTabs), replace = TRUE)
  
  # Subset for a single sample
  singleSampleSegTabs <- subset(segTabs, sample == "TCGA-A8-A07S-01A-11D-A036-01")
  
  # Read copy number data
  cn <- read_copynumber(singleSampleSegTabs,
                        seg_cols = c("chromosome", "start", "end", "segVal"),
                        genome_measure = "wg", complement = TRUE, add_loh = TRUE)
  
  # Tally with Wang method
  tally_wang <- sigminer::sig_tally(cn, method = "Wang", keep_only_matrix = FALSE)
  
  # Check that row names are not NULL
  expect_false(is.null(rownames(tally_wang$nmf_matrix)))
  
  # Check that row names contain the sample ID
  expect_equal(rownames(tally_wang$nmf_matrix), "TCGA-A8-A07S-01A-11D-A036-01")
  
  # For comparison, test Steele method as well
  tally_steele <- sigminer::sig_tally(cn, method = "S", keep_only_matrix = FALSE)
  
  # Both should have the same row names
  expect_equal(rownames(tally_wang$nmf_matrix), rownames(tally_steele$nmf_matrix))
})

test_that("sig_tally with Wang method preserves row names for multiple samples", {
  # Load test data with multiple samples
  load(system.file("extdata", "toy_copynumber.RData",
                   package = "sigminer", mustWork = TRUE))
  
  # Tally with Wang method
  cn_tally_W <- sig_tally(cn, method = "W", keep_only_matrix = FALSE)
  
  # Check that row names are not NULL
  expect_false(is.null(rownames(cn_tally_W$nmf_matrix)))
  
  # Check that we have row names for all samples
  expect_true(length(rownames(cn_tally_W$nmf_matrix)) > 0)
  
  # Check that row names are unique
  expect_equal(length(rownames(cn_tally_W$nmf_matrix)), 
               length(unique(rownames(cn_tally_W$nmf_matrix))))
})