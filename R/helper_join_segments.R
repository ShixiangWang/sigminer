
# data.table implementation 20% speed up ----------------------------------

helper_join_segments2 <- function(segTab) {
  segTab <- segTab[order(segTab$start)]
  final_orders <- c("chromosome", "start", "end", "segVal", "sample")
  cls_cols <- setdiff(colnames(segTab), c(final_orders, "segVal2"))
  segTab$segVal2 <- segTab$segVal
  segTab <- segTab[, .groupby_collapse(.SD, cls_cols),
    by = list(sample, chromosome)
  ]

  segTab$data.table <- NULL
  data.table::setcolorder(segTab, final_orders)
  return(segTab)
}

.groupby_collapse <- function(dt, cols = NULL) {
  if (length(cols) > 0) {
    dt <- dt[, collapse_top2bottom(.SD, cols = cols),
      by = data.table::rleid(segVal)
    ]
  } else {
    dt <- dt[, collapse_top2bottom(.SD),
      by = data.table::rleid(segVal)
    ]
  }
  dt
}

collapse_top2bottom <- function(dt, cols = NULL) {
  x <- data.table::data.table(
    start = dt[1]$start,
    end = dt[.N]$end,
    segVal = dt[1]$segVal2
  )

  if (!is.null(cols)) {
    x[, (cols) := lapply(dt[, cols, with = FALSE], function(x) {
      if (is.numeric(x)) {
        mean(x, na.rm = TRUE)
      } else {
        paste0(unique(na.omit(x)), collapse = ",")
      }
    })]
  }
  return(x)
}


# tidyverse implementation ------------------------------------------------


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

    cutpoint <- equal_index[which(diff(equal_index) > 1)] + 2L
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
          dplyr::summarise_at(
            res, dplyr::vars(-c("start", "end", "segVal")),
            ~ ifelse(is.numeric(.), mean(., na.rm = TRUE),
              paste0(unique(na.omit(.)), collapse = ",")
            )
          )
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
