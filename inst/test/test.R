library(maftools)
?read.maf
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read.maf(maf = laml.maf)
laml.tnm <- trinucleotideMatrix(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
                                prefix = 'chr', add = TRUE, useSyn = TRUE)
laml.sign <- extractSignatures(mat = laml.tnm, n = 3, pConstant = 0.001)
laml.sign$signatures

compareSignatures(laml.sign)
