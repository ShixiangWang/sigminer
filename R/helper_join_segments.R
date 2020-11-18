
# data.table implementation ~2x speed up ----------------------------------

helper_join_segments2 <- function(segTab, add_loh = FALSE, loh_min_frac = 0.3) {
  x <- segTab # This is just a alias, not a copy, modify inplace will save memory
  index <- seq(1, nrow(x))
  x$index <- index
  x$dval <- c(NA, diff(x$segVal))
  i_eq <- x[, list(index = .SD[-1]$index[.SD[-1]$dval == 0]),
    by = list(sample, chromosome)
  ]$index

  x$index <- NULL
  x$dval <- NULL

  if (length(i_eq) > 0) {
    ## Split joining segments into different groups
    i_join <- sort(union(i_eq - 1L, i_eq))
    # z <- split(i_eq, findInterval(i_eq, i_eq[diff(i_eq) > 1] + 2L)) ## the first element is not included
    # grp <- rep(names(z), sapply(z, length) + 1L)
    ## Use rle to improve performance
    z <- rle(findInterval(i_eq, i_eq[diff(i_eq) > 1] + 2L))
    z$lengths <- z$lengths + 1L
    z$values <- as.character(z$values)
    dt_join <- x[i_join][,
      .collapse_segTab(.SD,
        cols = setdiff(
          colnames(x),
          c("chromosome", "start", "end", "segVal", "sample")
        ),
        add_loh = add_loh,
        loh_min_frac = loh_min_frac
      ),
      by = list(.grp = inverse.rle(z))
    ]
    dt_join$.grp <- NULL

    x <- data.table::rbindlist(
      list(x[setdiff(index, i_join)], dt_join),
      fill = TRUE
    )
  }

  return(x)
}

.collapse_segTab <- function(dt, cols = NULL, add_loh = FALSE, loh_min_frac = 0.3) {
  x <- data.table::data.table(
    chromosome = dt[1]$chromosome,
    start = dt[1]$start,
    end = dt[.N]$end,
    segVal = dt[1]$segVal,
    sample = dt[1]$sample
  )

  if (length(cols) > 0) {
    x[, (cols) := lapply(dt[, cols, with = FALSE], function(x) {
      if (is.numeric(x)) {
        mean(x, na.rm = TRUE)
      } else if (is.logical(x)) {
        any(x, na.rm = TRUE)
      } else {
        paste0(unique(na.omit(x)), collapse = ",")
      }
    })]
  }

  if (add_loh) {
    # When LOH regions have same total copy number values as adjacent
    # regions, only label the segments harbor LOH with minimal length fraction
    # A logical column 'loh' should exist
    if (any(dt$loh)) {
      dt$.len <- dt$end - dt$start + 1
      x$.loh_frac <- sum(dt$.len[dt$loh]) / sum(dt$.len)
      # x$loh is already TRUE
      # reset it
      if (x$.loh_frac <= loh_min_frac) {
        x$loh <- FALSE
      }
    } else {
      x$.loh_frac <- NA_real_
    }
  }

  return(x)
}

# helper_join_segments2_old <- function(segTab) {
#   # segTab <- segTab[order(segTab$start)]
#   final_orders <- c("chromosome", "start", "end", "segVal", "sample")
#   cls_cols <- setdiff(colnames(segTab), c(final_orders, "segVal2"))
#   segTab$segVal2 <- segTab$segVal
#   segTab <- segTab[, .groupby_collapse(.SD, cls_cols),
#                    by = list(sample, chromosome)
#   ]
#
#   segTab$data.table <- NULL
#   data.table::setcolorder(segTab, final_orders)
#   return(segTab)
# }
#
# .groupby_collapse <- function(dt, cols = NULL) {
#   if (length(cols) > 0) {
#     dt <- dt[, .collapse_top2bottom(.SD, cols = cols),
#       by = data.table::rleid(segVal)
#     ]
#   } else {
#     dt <- dt[, .collapse_top2bottom(.SD),
#       by = data.table::rleid(segVal)
#     ]
#   }
#   dt
# }
#
# .collapse_top2bottom <- function(dt, cols = NULL) {
#   x <- data.table::data.table(
#     start = dt[1]$start,
#     end = dt[.N]$end,
#     segVal = dt[1]$segVal2
#   )
#
#   if (!is.null(cols)) {
#     x[, (cols) := lapply(dt[, cols, with = FALSE], function(x) {
#       if (is.numeric(x)) {
#         mean(x, na.rm = TRUE)
#       } else {
#         paste0(unique(na.omit(x)), collapse = ",")
#       }
#     })]
#   }
#   return(x)
# }


# tidyverse implementation ------------------------------------------------


# helper_join_segments <- function(segTab) {
#   final_orders <- c("chromosome", "start", "end", "segVal", "sample")
#   segTab <- segTab %>%
#     dplyr::as_tibble() %>%
#     dplyr::group_by(.data$sample, .data$chromosome) %>%
#     tidyr::nest() %>%
#     dplyr::summarise(data = purrr::map(.data$data, join_segments)) %>%
#     tidyr::unnest("data") %>%
#     dplyr::select(final_orders, dplyr::everything()) %>%
#     data.table::as.data.table()
#
#   segTab
# }
#
# join_segments <- function(df) {
#   df <- df %>% dplyr::arrange(.data$start)
#   equal_index <- which(diff(df$segVal) == 0)
#
#   to_join_index <- sort(union(equal_index, equal_index + 1))
#
#   if (length(to_join_index) > 0) {
#     not_join_index <- setdiff(1:nrow(df), to_join_index)
#
#     cutpoint <- equal_index[which(diff(equal_index) > 1)] + 2L
#     join_list <- split(
#       to_join_index,
#       findInterval(to_join_index, cutpoint)
#     )
#
#
#     join_df <- purrr::map_df(join_list, function(index, df) {
#       res <- df[index, ]
#       out <- dplyr::tibble(
#         start = res$start[1],
#         end = res$end[length(res$end)],
#         segVal = res$segVal[1]
#       )
#       if (ncol(res) == 3L) {
#         return(out)
#       } else {
#         dplyr::bind_cols(
#           out,
#           dplyr::summarise_at(
#             res, dplyr::vars(-c("start", "end", "segVal")),
#             ~ ifelse(is.numeric(.), mean(., na.rm = TRUE),
#               paste0(unique(na.omit(.)), collapse = ",")
#             )
#           )
#         )
#       }
#     }, df = df)
#
#     if (length(not_join_index) > 0) {
#       df <- dplyr::bind_rows(
#         join_df,
#         df[not_join_index, ]
#       ) %>%
#         dplyr::arrange(.data$start)
#
#       return(df)
#     } else {
#       return(join_df)
#     }
#   } else {
#     return(df)
#   }
# }
