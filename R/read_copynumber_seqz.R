#' Read Absolute Copy Number Profile from Sequenza Result Directory
#'
#' @param target_dir a directory path.
#' @param return_df if `TRUE`, return a `data.frame` directly, otherwise return a
#' [CopyNumber] object.
#' @param ... other parameters passing to [read_copynumber()].
#'
#' @return a `data.frame` or a `CopyNumber` object.
#' @export
read_copynumber_seqz <- function(target_dir, return_df = FALSE, ...) {
  stopifnot(length(target_dir) == 1L, dir.exists(target_dir))

  SAMPLE <- dir(target_dir, pattern = "_segments.txt")

  if (length(SAMPLE) < 1) {
    stop("Cannot find any sequenza result file, please check your path!", call. = FALSE)
  }

  SAMPLE <- sub("_segments.txt", "", SAMPLE)
  res <- purrr::map2_df(file.path(
    target_dir,
    paste(SAMPLE,
      "_segments.txt",
      sep = ""
    )
  ), SAMPLE, function(x, y) {
    message("Processing ", y)
    df <- data.table::fread(x, data.table = FALSE)
    df <- df %>%
      dplyr::select(c("chromosome", "start.pos", "end.pos", "CNt", "B")) %>%
      dplyr::mutate(sample = y)
    colnames(df) <- c("Chromosome", "Start.bp", "End.bp", "modal_cn", "minor_cn", "sample")
    df
  })

  if (return_df) {
    return(res)
  } else {
    return(read_copynumber(res, ...))
  }
}
