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
#' @param TD_size_cutoff a length-3 numeric vector used to specify the start, midpoint, end
#' segment size for determining tandem duplication size range, midpoint is used to split
#' TD into short TD and long TD. Default is 1Kb to 100Kb for short TD, 100Kb to 2Mb for long
#' TD.
#' @param TD_cn_cutoff a number defining the maximum copy number of TD,
#' default is `Inf`, i.e. no cutoff.
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
#' - TDP_pnas: tandem duplication phenotype score from `https://www.pnas.org/doi/10.1073/pnas.1520010113`,
#' the threshold `k` in reference is omitted.
#' \deqn{TDP = - \frac{\sum_{chr} |TD_{obs}-TD_{exp}|}{TD_{total}}}
#' where \eqn{TD_{total}} is the number of TD, \eqn{TD_{obs}} and
#' \eqn{TD_exp} are observed number of TD and expected number of TD for each chromosome.
#' - TDP: tandem duplication score used defined by our group work,
#' TD represents segment with copy number greater than 2.
#' \deqn{TD = \frac{TD_{total}}{\sum_{chr} |TD_{obs}-TD_{exp}|+1}}
#' - sTDP: TDP score for short TD.
#' - lTDP: TDP score for long TD.
#' - TDP_size : TDP region size (Mb).
#' - sTDP_size: sTDP region size (Mb).
#' - lTDP_size: lTDP region size(Mb).
#' - Chromoth_state: chromothripsis state score,
#' according to reference \doi{10.1016/j.cell.2013.02.023},
#' chromothripsis frequently leads to massive loss of segments on
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
#' d <- scoring(cn)
#' d
#'
#' d2 <- scoring(cn, TD_cn_cutoff = 4L)
#' d2
#' @testexamples
#' expect_s3_class(d, "data.frame")
#' expect_s3_class(d2, "data.frame")
scoring <- function(object, TD_size_cutoff = c(1e3, 1e5, 2e6),
                    TD_cn_cutoff = Inf) {
  stopifnot(
    inherits(object, "CopyNumber"), length(TD_size_cutoff) == 3,
    is.numeric(TD_size_cutoff), is.numeric(TD_cn_cutoff)
  )

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

  ## Get TD(P) scores
  ## TDP for tandem duplication phenotype from https://www.pnas.org/content/113/17/E2373
  ## shortTD and longTD self defined tandem duplication phenotype score
  dat_TD <- scoring_TD(data, TD_size_cutoff, TD_cn_cutoff)

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

  dat[order(dat$cnaBurden)]
}

scoring_TD <- function(data, TD_size_cutoff, TD_cn_cutoff) {
  autosomes <- paste0("chr", 1:22)
  ## Fill autosomes
  data_full <-
    data %>%
    dplyr::as_tibble() %>%
    dplyr::filter(.data$segVal <= TD_cn_cutoff) %>%
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
    chr_TD = sum(segVal > 2 & segLen >= TD_size_cutoff[1] & segLen <= TD_size_cutoff[3]),
    chr_sTD = sum(segVal > 2 & segLen >= TD_size_cutoff[1] & segLen <= TD_size_cutoff[2]),
    len_TD = sum(segLen[segVal > 2 & segLen >= TD_size_cutoff[1] & segLen <= TD_size_cutoff[3]]),
    len_sTD = sum(segLen[segVal > 2 & segLen >= TD_size_cutoff[1] & segLen <= TD_size_cutoff[2]])
  ), by = list(sample, chromosome)]
  data_chr$chr_lTD <- data_chr$chr_TD - data_chr$chr_sTD
  data_chr$len_lTD <- data_chr$len_TD - data_chr$len_sTD

  calc_TDP <- function(x) {
    sum(x) / (sum(abs(x - sum(x) / 22)) + 1)
  }

  calc_TDP_pnas <- function(x) {
    # https://www.pnas.org/content/113/17/E2373
    # k is omit here
    -sum(abs(x - sum(x) / 22)) / sum(x)
  }

  data_chr[, list(
    TD_count = sum(chr_TD),
    TDP = calc_TDP(chr_TD),
    sTDP = calc_TDP(chr_sTD),
    lTDP = calc_TDP(chr_lTD),
    TDP_size = sum(len_TD) / 1e6,
    sTDP_size = sum(len_sTD) / 1e6,
    lTDP_size = sum(len_lTD) / 1e6,
    TDP_pnas = calc_TDP_pnas(chr_TD)
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
      Chromoth_state = sum(count^2)
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


utils::globalVariables(
  c(
    "chr_TD", "chr_lTD", "chr_sTD", "segLen",
    "len_TD", "len_lTD", "len_sTD"
  )
)
