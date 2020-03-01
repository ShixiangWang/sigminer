#' Score Copy Number Profile
#'
#' Returns quantification of copy number profile and events including
#' tandem duplication and Chromothripisis etc.
#' Only copy number data from autosome is used here.
#' **Some of the quantification methods are rough,
#' you use at your risk**. You should do some extra work to check the
#' result scores.
#'
#' @param object a object of [CopyNumber].
#'
#' @return a `data.table` with following scores:
#' - cnaBurden: CNA burden representing the altered genomic fraction as previously reported.
#' - cnaLoad: CNA load representing the quantity of copy number alteration.
#' - MACN: mean altered copy number (MACN) reflecting the property of altered copy number segments,
#' calculated as
#' \deqn{MACN = \frac{\sum_{i} CN_i}{N_{cnv}}}
#' where \eqn{CN_i} is the copy number of altered segment \eqn{i}, \eqn{N_{cnv}} is
#' the number of CNV.
#' - weightedMACN: same as MACN but weighted with segment length.
#' \deqn{MACN_{weighted} = \frac{\sum_{i} (CN_i \times L_{i})}{ \sum_{i} L_{i} }}
#' where \eqn{L_{i}} is the length of altered copy number segment \eqn{i}.
#' - Ploidy: ploidy, the formula is same as `weightedMACN` but using all copy number segments instead of
#' altered copy number segments.
#' - TDP: tandem duplication phenotype score from <https://www.pnas.org/content/113/17/E2373>,
#' the threshold `k` in reference is omitted.
#' \deqn{TDP = - \frac{\sum_{chr} |TD_{obs}-TD_{exp}|}{TD_{total}}}
#' where \eqn{TD_{total}} is the number of TD, \eqn{TD_{obs}} and
#' \eqn{TD_exp} are observed number of TD and expected number of TD for each chromosome.
#' - TD: tandem duplication score, TD represents segment with copy number greater than 2.
#' \deqn{TD = \frac{TD_{total}}{\sum_{chr} |TD_{obs}-TD_{exp}|+1}}
#' - sTD: short tandem duplication score, same formula as `TD` but focus on segments with length less than or equal to 2Mb.
#' - lTD: long tandem duplication score, same formula as `TD` but focus on segments with length greater than 2Mb.
#' - Chromothripisis: according to reference <http://dx.doi.org/10.1016/j.cell.2013.02.023>,
#' Chromothripsis frequently leads to massive loss of segments on
#' the affected chromosome with segmental losses being interspersed with regions displaying
#' normal (disomic) copy-number (e.g., copy-number states oscillating between
#' copy-number = 1 and copy-number = 2), form tens to hundreds of locally clustered DNA rearrangements.
#' Most of methods use both SV and CNV to infer chromothripsis, here we roughly quantify it with
#' \deqn{\sum_{chr}{N_{OsCN}^2}}
#' where \eqn{N_{OsCN}} is the number of oscillating copy number pattern "2-1-2" for each chromosome.
#' @export
#'
#' @examples
#' # Load copy number object
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' scoring(cn)
scoring <- function(object) {
  stopifnot(inherits(object, "CopyNumber"))

  ## The data all from autosomes
  data <- object@data
  data <- data[data$chromosome %in% paste0("chr", 1:22)]
  data$segLen <- data$end - data$start + 1

  ## Get cnaLoad and cnaBurden
  summary <- object@summary.per.sample
  dat_cnaBurden <- summary[, c("sample", "cna_burden"), with = FALSE]
  colnames(dat_cnaBurden)[2] <- "cnaBurden"
  dat_cnaLoad <- data[, list(cnaLoad = round(log2(sum(segVal != 2) + 1), 2)), by = sample]
  ##


  ## Get Ploidy
  dat_Ploidy <- get_cn_ploidy(data)
  colnames(dat_Ploidy)[2] <- "Ploidy"
  ##

  ## Get MACN score
  ## MACN for mean altered copy number score
  ## also weighted by segment length
  dat_MACN <- scoring_MACN(data)
  ##

  ## Get TD scores
  ## TDP for tandem duplication phenotype from https://www.pnas.org/content/113/17/E2373
  ## shortTD and longTD self defined tandem duplication score, 2Mb as cutoff
  dat_TD <- scoring_TD(data)

  ## Get Chromothripisis score
  dat_Chromoth <- scoring_Chromoth(data)
  ##

  ## Combine data and return
  dat <- purrr::reduce(list(
    dat_cnaBurden,
    dat_cnaLoad,
    dat_MACN,
    dat_Ploidy,
    dat_TD,
    dat_Chromoth
  ), merge, by = "sample", all = TRUE)

  dat
}

scoring_TD <- function(data) {
  autosomes <- paste0("chr", 1:22)
  ## Fill autosomes
  data_full <-
    data %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(.data$sample) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      data = purrr::map(.data$data, function(x) {
        left_chrs <- setdiff(autosomes, unique(x$chromosome))
        n_left <- length(left_chrs)
        if (n_left > 0) {
          dplyr::bind_rows(
            x,
            dplyr::tibble(
              chromosome = left_chrs,
              start = rep(0, n_left),
              end = rep(0, n_left),
              segVal = rep(2L, n_left),
              segLen = rep(0, n_left)
            )
          )
        } else {
          x
        }
      })
    ) %>%
    tidyr::unnest("data") %>%
    data.table::as.data.table()

  data_chr <- data_full[, list(
    chr_TD = sum(segVal > 2),
    chr_sTD = sum(segVal > 2 & segLen <= 2e6)
  ), by = list(sample, chromosome)]
  data_chr$chr_lTD <- data_chr$chr_TD - data_chr$chr_sTD

  calcTD <- function(x) {
    sum(x) / (sum(abs(x - sum(x) / 22)) + 1)
  }

  calcTDP <- function(x) {
    # https://www.pnas.org/content/113/17/E2373
    # k is omit here
    -sum(abs(x - sum(x) / 22)) / sum(x)
  }

  ## Calculate scores
  # data_chr %>%
  #   dplyr::as_tibble() %>%
  #   dplyr::group_by(.data$sample) %>%
  #   dplyr::summarise(
  #     TDP = calcTDP(.data$chr_TD),
  #     TD = calcTD(.data$chr_TD),
  #     sTD = calcTD(.data$chr_sTD),
  #     lTD = calcTD(.data$chr_lTD)
  #   )

  data_chr[, list(
    TDP = calcTDP(chr_TD),
    TD = calcTD(chr_TD),
    sTD = calcTD(chr_sTD),
    lTD = calcTD(chr_lTD)
  ), by = list(sample)]
}

scoring_Chromoth <- function(data) {
  count <- function(x) {
    seqs_mid <- as.integer(x)
    seqs_left <- as.integer(dplyr::lag(x, default = 2))
    seqs_right <- as.integer(dplyr::lead(x, default = 2))
    seqs <- paste0(seqs_left, seqs_mid, seqs_right)
    sum(seqs == "212")
  }

  data %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(.data$sample, .data$chromosome) %>%
    dplyr::summarise(
      count = count(.data$segVal)
    ) %>%
    dplyr::summarise(
      Chromothripisis = sum(count^2)
    ) %>%
    data.table::as.data.table()
}


scoring_MACN <- function(data) {
  data %>%
    dplyr::as_tibble() %>%
    dplyr::filter(.data$segVal != 2) %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::summarise(
      MACN = sum(.data$segVal) / dplyr::n(),
      weightedMACN = sum(.data$segVal * .data$segLen) / sum(.data$segLen)
    ) %>%
    dplyr::mutate(
      MACN = round(.data$MACN, 2),
      weightedMACN = round(.data$weightedMACN, 2)
    ) %>%
    data.table::as.data.table()
}
