#' Obtain Signature Index for Cancer Types
#'
#' @param keyword keyword to search in the signature index database.
#' @param sig_type signature type.
#' @param seq_type sequencing type.
#' @param source data source.
#'
#' @return a `list`.
#' @export
#'
#' @examples
#' l1 <- get_sig_cancer_type_index()
#' l2 <- get_sig_cancer_type_index(sig_type = "SBS")
#' l3 <- get_sig_cancer_type_index(sig_type = "DBS", source = "PCAWG", seq_type = "WGS")
#' l4 <- get_sig_cancer_type_index(sig_type = "ID")
#' l5 <- get_sig_cancer_type_index(keyword = "breast")
#' l1
#' l2
#' l3
#' l4
#' l5
#' @testexamples
#' expect_is(l1, "list")
#' expect_is(l2, "list")
#' expect_is(l3, "list")
#' expect_is(l4, "list")
#' expect_null(l5)
get_sig_cancer_type_index <- function(
                                      sig_type = c("legacy", "SBS", "DBS", "ID"),
                                      seq_type = c("WGS", "WES"),
                                      source = c("PCAWG", "TCGA", "nonPCAWG"),
                                      keyword = NULL) {
  if (!is.null(keyword)) {
    df1 <- readRDS(system.file("extdata", "cosmic2_record_by_cancer.rds", package = "sigminer"))
    df2 <- readRDS(system.file("extdata", "signature_record_by_cancer.rds", package = "sigminer"))
    colnames(df1)[1] <- "name"
    colnames(df2)[1] <- "name"

    stopifnot(length(keyword) == 1L)

    keyword <- as.character(keyword)
    df1 <- subset(df1, grepl(keyword, df1$name, ignore.case = TRUE))
    df2 <- subset(df2, grepl(keyword, df2$name, ignore.case = TRUE))

    if (nrow(df1) > 0) {
      message("Info found in COSMIC v2 database:")
      colnames(df1) <- c("cancer type", "signature index")
      print(df1)
    } else {
      message("Nothing found in COSMIC v2 database.")
    }

    if (nrow(df2) > 0) {
      message("Info found in COSMIC v3 database:")
      colnames(df2)[1] <- "cancer type"
      df2 <- df2 %>%
        tidyr::pivot_longer(
          cols = c("SBS list", "DBS list", "ID list"),
          names_to = "signature type",
          values_to = "signature index"
        ) %>%
        as.data.frame()
      df2[["signature type"]] <- sub(" list", "", df2[["signature type"]])
      print(df2)
    } else {
      message("Nothing found in COSMIC v3 database.")
    }

    return(invisible(NULL))
  }

  sig_type <- match.arg(sig_type)
  seq_type <- match.arg(seq_type)
  data_source <- match.arg(source)

  if (sig_type == "legacy") {
    df <- readRDS(system.file("extdata", "cosmic2_record_by_cancer.rds", package = "sigminer"))

    if (requireNamespace("knitr", quietly = TRUE)) {
      print(knitr::kable(df))
    } else {
      print(df)
    }

    res <- list()
    for (i in df[["full name"]]) {
      res[[i]] <- split_seq(df[["signature type"]][df[["full name"]] == i])
    }
  } else {
    df <- readRDS(system.file("extdata", "signature_record_by_cancer.rds", package = "sigminer"))
    sel_cols <- c("cancer type", paste(sig_type, "list"))
    df <- df %>%
      dplyr::filter(.data$source == data_source & .data$scale == seq_type) %>%
      dplyr::select(sel_cols)

    if (nrow(df) < 1) {
      message("No records found, please try setting other options.")
      return(invisible(NULL))
    }

    if (requireNamespace("knitr", quietly = TRUE)) {
      print(knitr::kable(df))
    } else {
      print(df)
    }

    res <- list()
    for (i in df[["cancer type"]]) {
      res[[i]] <- split_seq(df[[2]][df[["cancer type"]] == i] %>% as.character())
    }
  }

  return(res)
}

split_seq <- function(x, sep = ",") {
  x <- unlist(strsplit(x, split = sep))
  gsub(" ", "", x)
}
