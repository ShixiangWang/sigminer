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

#---- get series
features = get_features(tcga_segTabs)
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
