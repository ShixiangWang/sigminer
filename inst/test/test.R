extdata_dir = system.file("extdata", package = "sigminer")
dir(extdata_dir)

cp = read_copynumber(extdata_dir, genome_build = "hg19", verbose = TRUE)
cp

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)
