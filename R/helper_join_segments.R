helper_join_segments <- function(segTab) {
  final_orders <- c("chromosome", "start", "end", "segVal", "sample")
  segTab <- segTab %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(.data$sample, .data$chromosome) %>%
    tidyr::nest() %>%
    dplyr::summarise(data = purrr::map(.data$data, join_segments)) %>%
    tidyr::unnest("data") %>%
    dplyr::select(final_orders, dplyr::everything()) %>%
    data.table::as.data.table()

  segTab
}

join_segments <- function(df) {
  df <- df %>% dplyr::arrange(.data$start)
  equal_index <- which(diff(df$segVal) == 0)

  to_join_index <- sort(union(equal_index, equal_index + 1))

  if (length(to_join_index) > 0) {
    not_join_index <- setdiff(1:nrow(df), to_join_index)

    cutpoint <- equal_index[which(diff(equal_index) > 1)] + 2
    join_list <- split(
      to_join_index,
      findInterval(to_join_index, cutpoint)
    )


    join_df <- purrr::map_df(join_list, function(index, df) {
      res <- df[index, ]
      out <- dplyr::tibble(
        start = res$start[1],
        end = res$end[length(res$end)],
        segVal = res$segVal[1]
      )
      if (ncol(res) == 3L) {
        return(out)
      } else {
        dplyr::bind_cols(
          out,
          dplyr::summarise_at(res, dplyr::vars(-c("start", "end", "segVal")),
                              ~ifelse(is.numeric(.), mean(., na.rm = TRUE),
                                      paste0(unique(na.omit(.)), collapse = ",")))
        )
      }
    }, df = df)

    if (length(not_join_index) > 0) {
      df <- dplyr::bind_rows(
        join_df,
        df[not_join_index, ]
      ) %>%
        dplyr::arrange(.data$start)

      return(df)
    } else {
      return(join_df)
    }
  } else {
    return(df)
  }
}
