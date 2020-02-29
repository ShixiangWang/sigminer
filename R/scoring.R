scoring = function(object) {
  stopifnot(inherits(data, "CopyNumber"))

  ## Get cnaLoad and cnaBurden
  summary = object@summary.per.sample
  summary$cnaLoad = log2(summary$n_of_cnv + 1)
  dat_CNA = summary[, c("sample", "cnaLoad", "cna_burden"), with = FALSE]
  colnames(dat_CNA)[3] = "cnaBurden"
  ##

  ## Get Ploidy
  dat_Ploidy = get_cn_ploidy(object)
  colnames(dat_Ploidy)[2] = "Ploidy"
  ##

  ## Get MACN score
  ## MACN for mean altered copy number score
  ## also weighted by segment length
  dat_MACN = scoring_MACN(object)
  ##

  cn_list <- get_cnlist(object, ignore_chrs = c("chrX", "chrY"))

  ## Get TD scores
  ## TDP for tandem duplication phenotype from https://www.pnas.org/content/113/17/E2373
  ## shortTD and longTD self defined tandem duplication score, 2Mb as cutoff
  dat_TD = scoring_TD(cn_list)

  ## Get Chromothripisis score
  dat_Chromoth = scoring_Chromoth(cn_list)
  ##

  ## Combine data and return
  ## dat

}

scoring_TD = function(object) {

  data = object@data
  data$segLen <- data$end - data$start + 1
  data_amp = data[data$segLen > 2]

  data_total = data_amp[, list(n_TD = .N,
                               n_sTD = sum(segLen <= 2e6)), by = sample]
  data_total$n_lTD = data_total$n_TD - data_total$n_sTD

  data_chr = data_amp[, list(chr_TD = .N,
                             chr_sTD = sum(segLen <= 2e6)), by = list(sample, chromosome)]
  data_chr$chr_lTD = data_chr$chr_TD - data_chr$chr_sTD

  ## Generate data for chr not exist


  data = merge(data_chr, data_total, by = "sample")

}

scoring_Chromoth = function(cn_list) {

}


scoring_MACN = function(object) {

  data = object@data
  data$segLen <- data$end - data$start + 1

  data %>%
    dplyr::as_tibble() %>%
    dplyr::filter(.data$segVal != 2) %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::summarise(
      MACN = sum(.data$segVal) / dplyr::n(),
      weightedMACN = sum(.data$segVal * .data$segLen) / sum(.data$segLen)
    ) %>%
    dplyr::mutate(
      MACN = round(.data$MACN, 2),
      weightedMACN = round(.data$weightedMACN, 2)
    ) %>%
    data.table::as.data.table()
}

