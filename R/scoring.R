scoring = function(object) {
  stopifnot(inherits(data, "CopyNumber"))

  summary = object@summary.per.sample
  summary$cnaLoad = log2(summary$n_of_cnv + 1)
  #cnaLoad =
  #cnaBurden =

  cn_list <- get_cnlist(object, ignore_chrs = c("chrX", "chrY"))


  ## TDP for tandem duplication phenotype from https://www.pnas.org/content/113/17/E2373
  ## shortTD and longTD self defined tandem duplication score, 2Mb as cutoff
  ## MACN for mean altered copy number score
  ## Ploidy
  ## Chromothripisis
  ## cnaLoad for number of cna
  ## cnaBurden for fraction of genome altered
  score_types <- c("TDP", "shortTP", "longTP",
                   "MACN", "Ploidy", "Chromothripisis",
                   "cnaLoad", cnaBurden)

}

scoring_TD = function(cn_list) {

}

scoring_MACN = function(cn_list) {

}

scoring_Ploidy = function(cn_list) {

}

scoring_Chromoth = function(cn_list) {

}



