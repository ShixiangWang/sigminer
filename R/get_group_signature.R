#' Find Enriched Signature in Each Sample Group
#'
#' @param group_df a `data.frame` containing columns 'sample' and 'group',
#' it can be obtained from [get_groups()].
#' @inheritParams get_sig_exposure
#'
#' @return a `data.frame`
#' @export
#'
#' @examples
#' \donttest{
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Extract copy number signatures
#' library(NMF)
#' sig <- sig_extract(cn_prepare$nmf_matrix, 2, nrun = 10)
#'
#' group_df = get_groups(sig, match_consensus = TRUE)
#' get_group_signature(group_df, sig)
#' }
get_group_signature = function(group_df, Signature) {
  exposure_df = get_sig_exposure(Signature = Signature, type = 'relative')
  merge_df = dplyr::left_join(group_df, exposure_df, by='sample')
  ztable = merge_df %>%
    tidyr::gather(key = 'Signature', value = 'Exposure', dplyr::starts_with("Sig")) %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::top_n(1, .data$Exposure) %>%
    dplyr::ungroup() %>%
    {table(.$group, .$Signature)}

  res = dplyr::tibble()
  for (i in seq_len(ncol(ztable))) {
    group = names(which.max(ztable[, i]))
    sig = colnames(ztable)[i]
    res = dplyr::bind_rows(res, dplyr::tibble(group=group, enrich_sig=sig))
  }

  merge_df %>%
    dplyr::select(c('sample', "group")) %>%
    dplyr::left_join(res, by = c('group')) %>%
    dplyr::as_tibble()
}
