extdata_dir = system.file("extdata", package = "sigminer")
dir(extdata_dir)

#---------------- Reading ----------------------------#

#------------ Test reading as CopyNumber object
cp = read_copynumber(extdata_dir, pattern = "txt", genome_build = "hg19", verbose = TRUE)
cp = read_copynumber(extdata_dir, pattern = "txt", genome_build = "hg19")
cp

#------------ Test reading as MAF object
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)

#------------ Test reading as GenomicVariation object
vt = read_variation(cp,  laml)
vt

#-------------- NMF prepare for Copy Number ------------#
example_data = system.file("extdata", "example_cn_list.RData", package = "sigminer")
load(example_data)

tcga_segTabs2 = lapply(tcga_segTabs,
                       function(x) {
                         data.table::as.data.table(x)
                         x$chromosome = paste0("chr", x$chromosome)
                         x = subset(x, chromosome %in% paste0("chr", 1:22))
                         data.table::as.data.table(x)})

#---- get series
features = get_features(tcga_segTabs2)
components = get_components(features)
inp_matrix = get_matrix(features, components)

#----
cn_list = get_cnlist(cp)

cn_features = get_features(cn_list, cores = 1)
cn_components = get_components(cn_features, nrep = 20)
cn_matrix = get_matrix(cn_features, cn_components, cores = 1)

cn_features = get_features(cn_list, cores = 2)
cn_matrix2 = get_matrix(cn_features, cn_components, cores = 2)

#--- prepare series
cn_matrix3 = prepare_copynumber(cp)
cn_matrix4 = prepare_copynumber(cp, reference_components = NULL)

#-------------- NMF prepare for MAF ------------#
laml.tnm <- prepare_maf(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
                        prefix = 'chr', add = TRUE, useSyn = TRUE)
laml.tnm2 = maftools::trinucleotideMatrix(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
                                          prefix = 'chr', add = TRUE, useSyn = TRUE)

#-------------- Signature Analysis Series --------------#

#- Prepare CopyNumber
nmf_matrix.cp = sig_prepare(cp)
#- Prepare MAF
nmf_matrix.maf = sig_prepare(laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
                             prefix = 'chr', add = TRUE, useSyn = TRUE)

#- Signature best rank estimation
sig_estimate(laml.tnm$nmf_matrix)
sig_estimate(inp_matrix, verbose = TRUE)

#------------- Test Maftools ----------------
require(maftools)

laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') #path to TCGA LAML MAF file
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') # clinical information containing survival information and histology. This is optional

laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
laml.tnm = trinucleotideMatrix(maf = laml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
library('NMF')
laml.sign = extractSignatures(mat = laml.tnm, nTry = 6, plotBestFitRes = FALSE)


#########################################
getBPnum2 <- function(abs_profiles, chrlen)
{
  out <- c()
  samps <- names(abs_profiles)
  for (i in samps)
  {
    message("Processing samps ", i)
    segTab <- abs_profiles[[i]]

    chrs <- unique(segTab$chromosome)

    allBPnum <- c()
    for (c in chrs)
    {
      currseg <- segTab[chromosome == c,]
      message("Processing chrs ", c)
      intervals <-
        seq(1, chrlen[chrlen[, 1] == c, 2] + 10000000, 10000000)
      res <-
        hist(currseg$end[-nrow(currseg)],
             breaks = intervals,
             plot = FALSE)$counts
      allBPnum <- c(allBPnum, res)
    }
    out <-
      rbind(out, cbind(ID = rep(i, length(allBPnum)), value = allBPnum))
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}

getBPnum2(tcga_segTabs2, chrlen = sigminer::chromsize.hg19)
