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

tcga_segTabs2 = lapply(tcga_segTabs, data.table::as.data.table)

features = get_features(tcga_segTabs)
components = get_components(features)
inp_matrix = get_matrix(features, components)

