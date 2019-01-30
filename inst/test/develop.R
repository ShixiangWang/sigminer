##########################################
### Store temp thoughts and functions ####
##########################################

load(system.file("extdata", "example_cn_list.RData",
                 package = "sigminer", mustWork = TRUE))

segTabs = data.table::rbindlist(tcga_segTabs, idcol = "sample")
samps = segTabs[, .N, by=sample][order(N)]
segTabs = segTabs[sample %in% samps$sample[6:15]]

cn = read_copynumber(segTabs,
                     seg_cols = c("chromosome", "start", "end", "segVal"),
                     genome_build = "hg19", complement = FALSE, verbose = TRUE)
read_copynumber(segTabs,
                     seg_cols = c("chromosome", "start", "end", "segVal"),
                     genome_build = "hg19", genome_measure = "wg", verbose = T)


cn_prepare =  sig_prepare(cn)

library(NMF)
sig_estimate(cn_prepare$nmf_matrix, verbose = TRUE)
res = sig_extract(cn_prepare$nmf_matrix, 2, mode = "copynumber", nrun = 1, cores = 1)

load(system.file("extdata", "toy_copynumber_signature.RData",
                             package = "sigminer", mustWork = TRUE))

subtypes = sig_assign_samples(res$nmfObj, type = "samples")
set.seed(1234)
subtypes$new_group = sample(c("1", "2","3", "4"), size = nrow(subtypes), replace = TRUE)
subtypes.sum = sig_summarize_subtypes(subtypes[, -1], col_subtype = "nmf_subtypes",
                                     cols_to_summary = colnames(subtypes[, -1])[c(-1,-2)],
                                     type = c("co", "ca"), verbose = TRUE)


draw_sig_profile(res$nmfObj)
draw_sig_activity(res$nmfObj)

sig_activity = sig_get_activity(res$nmfObj)
sig_cor = sig_get_correlation(sig_activity)
