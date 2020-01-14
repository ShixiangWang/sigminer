# Find Enriched Signature in Each Sample Group Based on Its Dominant Signature (Weighted by Exposure)
find_enriched_signature <- function(group_df, Signature) {

  # Handle hyper mutant samples
  hyper_index <- grepl("_\\[hyper\\]_", group_df$sample)
  if (sum(hyper_index) > 0) {
    sample.hyper <- sapply(
      group_df$sample,
      function(x) strsplit(x, "_\\[hyper\\]_")[[1]][[1]]
    )

    df <- dplyr::tibble(
      sample = names(sample.hyper),
      ID = as.character(sample.hyper)
    )

    group_df <- dplyr::left_join(group_df, df, by = "sample") %>%
      dplyr::select(-.data$sample) %>%
      dplyr::rename(sample = .data$ID) %>%
      dplyr::group_by(.data$sample) %>%
      dplyr::summarise(
        group = table(.data$group) %>% sort(decreasing = TRUE) %>% names() %>% .[[1]],
        silhouette_width = mean(.data$silhouette_width)
      )
  }

  exposure_df <- get_sig_exposure(Signature = Signature, type = "relative")
  merge_df <- dplyr::left_join(group_df, exposure_df, by = "sample")
  ztable <- merge_df %>%
    tidyr::gather(key = "Signature", value = "Exposure", dplyr::starts_with("Sig")) %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::top_n(1, .data$Exposure) %>%
    dplyr::ungroup() %>%
    {
      table(.$group, .$Signature)
    }

  message("=> Generating a table of group and dominant signature:")
  print(ztable)
  message("=> Assigning a group to a signature with the maxium fraction (stored in 'map_table' attr)...")

  res <- dplyr::tibble()
  for (i in seq_len(nrow(ztable))) {
    sig <- names(which.max(ztable[i, ]))
    res <- dplyr::bind_rows(res, dplyr::tibble(group = rownames(ztable)[i], enrich_sig = sig))
  }

  df <- merge_df %>%
    dplyr::select(-dplyr::starts_with("Sig")) %>%
    dplyr::left_join(res, by = c("group"))

  list(data = df, table = ztable)
}
