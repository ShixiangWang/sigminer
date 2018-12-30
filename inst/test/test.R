extdata_dir = system.file("extdata", package = "sigminer")
dir(extdata_dir)

#---------------- Reading ----------------------------#

#------------ Test reading as CopyNumber object
cp = read_copynumber(extdata_dir, genome_build = "hg19", verbose = TRUE)
cp

#------------ Test reading as MAF object
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)

#------------ Test reading as GenomicVariation object
vt = read_variation(cp,  laml)
vt

