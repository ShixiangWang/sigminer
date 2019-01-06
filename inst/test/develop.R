##########################################
### Store temp thoughts and functions ####
##########################################

load(system.file("extdata", "example_cn_list.RData",
                 package = "sigminer", mustWork = TRUE))

segTabs = data.table::rbindlist(tcga_segTabs, idcol = "sample")

cn = read_copynumber(segTabs,
                     seg_cols = c("chromosome", "start", "end", "segVal"),
                     genome_build = "hg19")


cn_prepare =  sig_prepare(cn, cores = 2)

res = sig_extract(cn_prepare$nmf_matrix, 3, mode = "copynumber", nrun = 5, cores = 2)

