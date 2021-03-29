#' Get Shannon Diversity Index for Signatures
#'
#' \deqn{H = - \sum_{i=1}^n{p_i ln(p_i)}}
#' where `n` is the number
#' of signatures identified in the signature with exposure > `cutoff`,
#' and `pi` is the normalized exposure of the ith signature with
#' exposure > `cutoff`. Exposures of signatures were normalized to
#' sum to `1`.
#'
#' @param rel_expo a `data.frame` with numeric columns indicating
#' **relative** signature exposures for each sample. Typically
#' this data can be obtained from [get_sig_exposure()].
#' @param cutoff a relative exposure cutoff for filtering signatures,
#' default is `0.1%`.
#'
#' @return a `data.frame`
#' @references
#' Steele, Christopher D., et al. "Undifferentiated sarcomas develop through distinct evolutionary pathways." Cancer Cell 35.3 (2019): 441-456.
#' @export
#'
#' @examples
#' # Load mutational signature
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Get signature exposure
#' rel_expo <- get_sig_exposure(sig2, type = "relative")
#' rel_expo
#' diversity_index <- get_shannon_diversity_index(rel_expo)
#' diversity_index
#' @testexamples
#' expect_is(rel_expo, "data.frame")
#' expect_is(diversity_index, "data.frame")
get_shannon_diversity_index <- function(rel_expo, cutoff = 0.001) {
  stopifnot(is.data.frame(rel_expo))
  dt <- rel_expo %>%
    dplyr::mutate_if(
      is.numeric,
      ~ ifelse(. > cutoff, -. * log(.), 0)
    ) %>%
    as.data.frame()
  num_idx <- sapply(dt, is.numeric)
  dplyr::bind_cols(
    dt[, !num_idx, drop = FALSE],
    data.frame(
      diversity_index = rowSums(dt[, num_idx, drop = FALSE])
    )
  ) %>%
    data.table::as.data.table()
}
