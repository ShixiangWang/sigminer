test_that("read_maf works", {
  laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
  laml <- read_maf(maf = laml.maf, verbose = FALSE)
  expect_s4_class(laml, "MAF")
})

test_that("read_copynumber works", {
  # Load toy dataset of absolute copynumber profile
  load(system.file("extdata", "toy_segTab.RData",
                   package = "sigminer", mustWork = TRUE
  ))
  cn <- read_copynumber(segTabs,
                        seg_cols = c("chromosome", "start", "end", "segVal"),
                        genome_build = "hg19", complement = FALSE, verbose = FALSE
  )
  expect_s4_class(cn, "CopyNumber")
  cn <- read_copynumber(segTabs,
                        seg_cols = c("chromosome", "start", "end", "segVal"),
                        genome_build = "hg19", complement = TRUE, verbose = FALSE
  )
  expect_s4_class(cn, "CopyNumber")
})
