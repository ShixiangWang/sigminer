#' Get Sample Groups from Signature Decomposition Information
#'
#' One of key results from signature analysis is to cluster samples into different
#' groups. This function takes `Signature` object as input
#' and return the membership in each cluster.
#'
#' Users may find there are bigger differences between using method 'samples' and 'exposure' but
#' they use a similar idear to find dominant signature, here goes the reason:
#'
#' Method 'samples' using data directly from NMF decomposition, this means the two matrix
#' `W` (basis matrix or signature matrix) and `H` (coefficient matrix or exposure matrix) are
#' the results of NMF. For method 'exposure', it uses the signature exposure loading matrix.
#' In this situation, each signture represents a number of mutations (alterations)
#  assigned to each signature, thus gives a more biologically meaningful result. More
#' about implementation please see source code of [sig_extract()] function.
#'
#' @param Signature a `Signature` object obtained either from [sig_extract] or [sig_auto_extract].
#' @param method grouping method, more see details, could be one of the following:
#' - 'consensus' - returns the cluster membership based on the hierarchical clustering of the consensus matrix,
#' it can only be used for the result obtained by [sig_extract()] with multiple runs using **NMF** package.
#' - 'k-means' -  returns the clusters by k-means.
#' - 'exposure' - assigns a sample into a group whose signature exposure
#' is dominant.
#' - 'samples' - returns the cluster membership based on the contribution of signature to each sample,
#' it can only be used for the result obtained by [sig_extract()] using **NMF** package.
#' @param n_cluster only used when the `method` is 'k-means'.
#' @param match_consensus only used when the `method` is 'consensus'.
#' If `TRUE`, the result will match order as shown in consensus map.
#' @return a `data.table` object
#' @import NMF
#' @importFrom stats kmeans
#' @export
#' @examples
#' \donttest{
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_tally_W.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Extract copy number signatures
#' library(NMF)
#' sig <- sig_extract(cn_tally_W$nmf_matrix, 2,
#'   nrun = 10,
#'   pConstant = 1e-13
#' )
#'
#' # Methods 'consensus' and 'samples' are from NMF::predict()
#' get_groups(sig, method = "consensus", match_consensus = TRUE)
#' get_groups(sig, method = "samples")
#'
#' # Use k-means clustering
#' get_groups(sig, method = "k-means")
#' }
#' @seealso [NMF::predict()]
get_groups <- function(Signature,
                       method = c("consensus", "k-means", "exposure", "samples"),
                       n_cluster = NULL,
                       match_consensus = TRUE) {
  stopifnot(inherits(Signature, "Signature"), is.null(n_cluster) | n_cluster > 1)
  method <- match.arg(method)

  if (method == "consensus") {
    message("=> Obtaining clusters from the hierarchical clustering of the consensus matrix...")
    if (!"nmf_obj" %in% names(Signature$Raw)) {
      stop("Input Signature object does not contain NMF object, please select other methods")
    }
    nmfObj <- Signature$Raw$nmf_obj
    predict.consensus <- predict(nmfObj, what = "consensus")
    silhouette.consensus <- silhouette(nmfObj, what = "consensus")

    data <- data.frame(
      sample = sampleNames(nmfObj),
      nmf_subgroup = predict.consensus,
      sil_width = signif(silhouette.consensus[, "sil_width"], 3),
      stringsAsFactors = FALSE
    )
    colnames(data) <- c("sample", "group", "silhouette_width")
    # If we want to display as we see in consensusmap, we just need to reoder everything.
    # Now re-order data to match consensusmap sample order
    if (match_consensus) {
      sample.order <- attributes(predict.consensus)$iOrd
      data <- data[sample.order, ]
    }

    data$group <- as.character(data$group)
    message("=> Finding the dominant signature of each group...")
    data <- find_enriched_signature(data, Signature)
    ztable <- data$table
    data <- data$data
  } else if (method == "samples") {
    message("=> Obtaining clusters by the contribution of signature to each sample...")

    if (!"nmf_obj" %in% names(Signature$Raw)) {
      stop("Input Signature object does not contain NMF object, please select other methods")
    }
    nmfObj <- Signature$Raw$nmf_obj
    predict.samples <- predict(nmfObj, what = "samples", prob = T)
    silhouette.samples <- silhouette(nmfObj, what = "samples")
    data <- data.frame(
      sample = names(predict.samples$predict),
      group = predict.samples$predict,
      silhouette_width = signif(silhouette.samples[, "sil_width"], 3),
      prob = signif(predict.samples$prob, 3),
      stringsAsFactors = FALSE
    )

    data$group <- as.character(data$group)
    message("=> Finding the dominant signature of each group...")
    data <- find_enriched_signature(data, Signature)
    ztable <- data$table
    data <- data$data
  } else if (method == "exposure") {
    message("=> Creating clusters by the dominant signature (fraction is returned as weight)...")
    expo_df <- get_sig_exposure(Signature, type = "relative")
    data <- expo_df %>%
      tidyr::gather(key = "Signature", value = "Exposure", dplyr::starts_with("Sig")) %>%
      dplyr::group_by(.data$sample) %>%
      dplyr::top_n(1, .data$Exposure) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        enrich_sig = .data$Signature
      ) %>%
      dplyr::rename(
        group = .data$Signature,
        weight = .data$Exposure
      ) %>%
      dplyr::mutate(group = as.integer(sub("Sig", "", .data$group))) %>%
      dplyr::arrange(.data$group)

    data$group <- as.character(as.integer(factor(data$group)))
    ztable <- table(data$group, data$enrich_sig)
  } else if (method == "k-means") {
    set.seed(seed = 1024)
    expo_df <- get_sig_exposure(Signature, type = "relative")
    contrib <- expo_df %>%
      as.data.frame() %>%
      tibble::column_to_rownames("sample")
    n_cluster <- ifelse(is.null(n_cluster), ncol(contrib), n_cluster)
    message("=> Running k-means with ", n_cluster, " clusters...")
    contrib.km <- kmeans(x = contrib, centers = n_cluster)
    message("=> Generating a table of group and signature contribution (stored in 'map_table' attr):")
    ztable <- contrib.km$centers
    print(ztable)
    message("=> Assigning a group to a signature with the maximum fraction...")
    cluster_df <- as.data.frame(apply(t(ztable), 2, function(x) which(x == max(x))))
    colnames(cluster_df)[1] <- "enrich_sig"
    cluster_df$enrich_sig <- colnames(contrib)[cluster_df$enrich_sig]
    data.table::setDT(x = cluster_df, keep.rownames = TRUE)
    colnames(cluster_df)[1] <- "group"
    data <- as.data.frame(contrib.km$cluster)
    colnames(data)[1] <- "group"
    data.table::setDT(data, keep.rownames = TRUE)
    colnames(data)[1] <- "sample"
    data$group <- as.character(data$group)
    data <- merge(data, cluster_df, by = "group")
    # Set a default value for now
    # data$weight <- 1L
    # data.table::setcolorder(data, neworder = c("sample", "group", "weight", "enrich_sig"))
    data.table::setcolorder(data, neworder = c("sample", "group", "enrich_sig"))
  }

  data <- data.table::as.data.table(data)
  if (!match_consensus) {
    data <- data[order(as.integer(data$group))]
  }
  message("=> Summarizing...")
  sum_tb <- table(data$group)
  map_dt <- unique(data[, c("group", "enrich_sig"), with = FALSE])
  map_dic <- map_dt$enrich_sig
  names(map_dic) <- map_dt$group
  message(paste(paste0(
    "\tgroup #", names(sum_tb), ": ",
    sum_tb, " samples with ",
    map_dic[names(sum_tb)], " enriched."
  ), collapse = "\n"))

  if (exists("ztable")) {
    attr(data, "map_table") <- ztable
  }

  return(data)
}

utils::globalVariables(
  c("silhouette")
)
