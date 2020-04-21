#' Performs Strand Bias Enrichment Analysis for a Given Sample-by-Component Matrix
#'
#' See [sig_tally] for examples.
#'
#' @param mat a sample-by-component matrix from [sig_tally] with strand bias labels "T:" and "B:".
#'
#' @return a `data.table` sorted by `p_value`.
#' @export
enrich_component_strand_bias <- function(mat) {
  stopifnot(is.matrix(mat))

  ## mat: sample-by-component matrix with T: and U: labels
  T_cols <- startsWith(colnames(mat), "T:")
  U_cols <- startsWith(colnames(mat), "U:")

  if (!any(T_cols) & !any(U_cols)) {
    stop("No transcription labels found!")
  }

  types <- colnames(mat)[T_cols]
  types <- sub("T:", "", types)
  mat <- mat[, T_cols | U_cols]

  dt <- dplyr::bind_cols(
    dplyr::tibble(sample = rownames(mat)),
    dplyr::as_tibble(mat)
  )

  f4row <- function(df) {
    df <- df %>%
      tidyr::pivot_longer(dplyr::contains(":"),
        names_sep = ":",
        names_to = c("strand", "component"),
        values_to = "count"
      ) %>%
      dplyr::mutate(
        strand = ifelse(.data$strand == "T", "Trans", "UnTrans")
      ) %>%
      tidyr::pivot_wider(
        id_cols = "component",
        names_from = "strand",
        values_from = "count"
      ) %>%
      dplyr::mutate(
        Trans_Total = sum(.data$Trans),
        UnTrans_Total = sum(.data$UnTrans),
        Enrichment = .data$Trans / .data$UnTrans
      )


    df$p_value <- apply(df, 1, function(x) {
      x <- as.numeric(x[2:5])
      # The order is
      # U_Total-U, T_Total-T, U, T
      fisher.test(
        matrix(
          c(
            x[4] - x[2], x[3] - x[1],
            x[2], x[1]
          ),
          nrow = 2
        )
      )$p.value
    })
    df$fdr <- p.adjust(df$p_value, method = "fdr")

    return(df)
  }

  res <- dt %>%
    dplyr::group_by(.data$sample) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      data = list(
        purrr::map_df(
          .data$data,
          f4row
        )
      )
    ) %>%
    tidyr::unnest("data") %>%
    ## Should adjust all samples for a type instead
    ## of adjust all types for a sample??
    # dplyr::group_by(.data$component) %>%
    # tidyr::nest() %>%
    # dplyr::mutate(
    #   data = list(
    #     purrr::map_df(
    #       .data$data,
    #       ~dplyr::mutate(.x, fdr = p.adjust(.data$p_value, method = "fdr"))
    #     )
    #   )
    # ) %>%
    # tidyr::unnest("data") %>%
    dplyr::arrange(.data$p_value, .data$fdr) %>%
    data.table::as.data.table()

  colnames(res)[3:6] <- c("N_T", "N_U", "Total_T", "Total_U")

  return(res)
}
