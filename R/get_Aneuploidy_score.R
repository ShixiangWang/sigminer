#' Get Aneuploidy Score from Copy Number Profile
#'
#' This implements a Cohen-Sharir method (see reference) like "Aneuploidy Score" computation.
#' You can read the source code to see how it works. Basically, it follows
#' the logic of Cohen-Sharir method but with some difference in detail implementation.
#' Their results should be counterpart, but with no data validation for now.
#' **Please raise an issue if you find problem/bugs in this function**.
#'
#' @inheritParams read_copynumber
#' @param data a CopyNumber object or a `data.frame` containing at least
#' 'chromosome', 'start', 'end', 'segVal', 'sample' these columns.
#' @param ploidy_df default is `NULL`, compute ploidy by segment-size weighted copy number
#' aross autosome, see [get_cn_ploidy]. You can also provide a `data.frame` with 'sample'
#' and 'ploidy' columns.
#' @param rm_black_arms if `TRUE`, remove short arms of chr13/14/15/21/22 from calculation
#' as documented in reference #3.
#' @references
#' - Cohen-Sharir, Y., McFarland, J. M., Abdusamad, M., Marquis, C., Bernhard, S. V., Kazachkova, M., ... & Ben-David, U. (2021). Aneuploidy renders cancer cells vulnerable to mitotic checkpoint inhibition. Nature, 1-6.
#' - Logic reference: <https://github.com/quevedor2/aneuploidy_score/>.
#' - Taylor, Alison M., et al. "Genomic and functional approaches to understanding cancer aneuploidy." Cancer cell 33.4 (2018): 676-689.
#'
#' @return A `data.frame`
#' @export
#'
#' @examples
#' # Load copy number object
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' df <- get_Aneuploidy_score(cn)
#' df
#'
#' df2 <- get_Aneuploidy_score(cn@data)
#' df2
#'
#' df3 <- get_Aneuploidy_score(cn@data,
#'   ploidy_df = get_cn_ploidy(cn@data)
#' )
#' df3
#' @testexamples
#' expect_equal(nrow(df), 10L)
#' expect_equal(nrow(df2), 10L)
#' expect_equal(nrow(df3), 10L)
get_Aneuploidy_score <- function(data, ploidy_df = NULL, genome_build = "hg19", rm_black_arms = FALSE) {
  stopifnot(is.data.frame(data) | inherits(data, "CopyNumber"))
  if (is.data.frame(data)) {
    nc_cols <- c("chromosome", "start", "end", "segVal", "sample")
    if (!all(nc_cols %in% colnames(data))) {
      stop("Invalid input, it must contain columns: ", paste(nc_cols, collapse = " "))
    }
    annot <- get_LengthFraction(data,
      genome_build = genome_build,
      seg_cols = c("chromosome", "start", "end", "segVal")
    )
  } else {
    if (nrow(data@annotation) == 0L) {
      annot <- get_LengthFraction(data@data,
        genome_build = data@genome_build,
        seg_cols = c("chromosome", "start", "end", "segVal")
      )
    } else {
      annot <- data@annotation
    }
    data <- data@data
  }

  if (is.null(ploidy_df)) {
    ploidy_df <- get_cn_ploidy(data)
  } else {
    ploidy_df <- data.table::as.data.table(ploidy_df)
  }

  annot <- merge(annot, ploidy_df, by = "sample", all.x = TRUE)
  annot$type <- sub("[0-9]+", "", annot$location)

  # For segments spans the centromere, remove it if length < 50%, i.e. 1 here.
  # We assign average weights to all arms.
  df <- annot[!(type == "pq" & round(fraction, 2) < 1) & (chromosome %in% paste0("chr", 1:22)), ] %>%
    dplyr::mutate(
      fraction = ifelse(type == "pq", .data$fraction / 2, .data$fraction),
      type = ifelse(type == "pq", "p,q", type)
    ) %>%
    tidyr::separate_rows(type, sep = ",") %>%
    data.table::as.data.table()

  if (rm_black_arms) {
    # Remove short arms of chr13/14/15/21/22
    df <- df[!(chromosome %in% c("chr13", "chr14", "chr15", "chr21", "chr22") & type %in% "p"), ]
  }

  df <- df[, .(flag = as.integer(round(sum(segVal * fraction) - ploidy[1]))),
    # overage = round(sum(fraction), 3)),
    by = .(sample, chromosome, type)
  ]

  df %>%
    dplyr::mutate(
      flag = dplyr::case_when(
        .data$flag > 1 ~ 1L,
        .data$flag < -1 ~ -1L,
        TRUE ~ .data$flag
      ),
      cytoband = paste(.data$chromosome, .data$type, sep = "-")
    ) %>%
    dplyr::select(-c("chromosome", "type")) %>%
    tidyr::pivot_wider(id_cols = "sample", names_from = "cytoband", values_from = "flag", values_fill = 0L) %>%
    dplyr::mutate(AScore = rowSums(abs(.[, -1]))) %>%
    dplyr::select(c("sample", "AScore", dplyr::everything())) %>%
    data.table::as.data.table()
}

utils::globalVariables(c("ploidy"))
