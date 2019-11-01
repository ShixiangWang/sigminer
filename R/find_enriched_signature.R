# Find Enriched Signature in Each Sample Group Based on Its Dominant Signature (Weighted by Exposure)
find_enriched_signature = function(group_df, Signature) {
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
    dplyr::select(-dplyr::starts_with("Sig")) %>%
    dplyr::left_join(res, by = c('group'))
}
