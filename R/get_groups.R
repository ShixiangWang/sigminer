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
#' Now it can be used to relative exposure result in `data.table` format from [sig_fit].
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
#'   nrun = 10
#' )
#'
#' # Methods 'consensus' and 'samples' are from NMF::predict()
#' g1 <- get_groups(sig, method = "consensus", match_consensus = TRUE)
#' g1
#' g2 <- get_groups(sig, method = "samples")
#' g2
#'
#' # Use k-means clustering
#' g3 <- get_groups(sig, method = "k-means")
#' g3
#' }
#' @testexamples
#' expect_is(g1, "data.frame")
#' expect_is(g2, "data.frame")
#' expect_is(g3, "data.frame")
#' @seealso [NMF::predict()], [show_groups].
get_groups <- function(Signature,
                       method = c("consensus", "k-means", "exposure", "samples"),
                       n_cluster = NULL,
                       match_consensus = TRUE) {
  fit_flag <- data.table::is.data.table(Signature)
  stopifnot(
    inherits(Signature, "Signature") | fit_flag,
    is.null(n_cluster) | n_cluster > 1
  )
  method <- match.arg(method)

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  if (fit_flag) {
    send_success("A 'data.table' detected.")

    if (method %in% c("consensus", "samples")) {
      send_stop("Method 'consensus' and 'samples' cannot be applied to input data in {.code data.table} format.\n Choose others please.")
    }
    send_success("Method checked.")

    if (purrr::map_lgl(Signature, ~ ifelse(is.numeric(.), any(. > 1), FALSE)) %>% any()) {
      send_stop("When input is {.code data.table} (from sig_fit), a relative exposure result is valid.")
    }
    send_success("Exposure should be relative checked.")
  } else {
    send_success("'Signature' object detected.")
  }

  if (method == "consensus") {
    send_info("Obtaining clusters from the hierarchical clustering of the consensus matrix...")
    if (!"nmf_obj" %in% names(Signature$Raw)) {
      send_stop("Input Signature object does not contain NMF object, please select other methods.")
    }
    nmfObj <- Signature$Raw$nmf_obj
    predict.consensus <- predict(nmfObj, what = "consensus")
    sil <- "NMF" %:::% "silhouette.NMF"
    silhouette.consensus <- sil(nmfObj, what = "consensus")

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
    send_info("Finding the dominant signature of each group...")
    data <- find_enriched_signature(data, Signature)
    ztable <- data$table
    data <- data$data
  } else if (method == "samples") {
    send_info("Obtaining clusters by the contribution of signature to each sample...")

    if (!"nmf_obj" %in% names(Signature$Raw)) {
      send_stop("Input Signature object does not contain NMF object, please select other methods.")
    }
    nmfObj <- Signature$Raw$nmf_obj
    predict.samples <- predict(nmfObj, what = "samples", prob = T)
    silhouette.samples <- cluster::silhouette(nmfObj, what = "samples")
    data <- data.frame(
      sample = names(predict.samples$predict),
      group = predict.samples$predict,
      silhouette_width = signif(silhouette.samples[, "sil_width"], 3),
      prob = signif(predict.samples$prob, 3),
      stringsAsFactors = FALSE
    )

    data$group <- as.character(data$group)
    send_info("Finding the dominant signature of each group...")
    data <- find_enriched_signature(data, Signature)
    ztable <- data$table
    data <- data$data
  } else if (method == "exposure") {
    send_info("Creating clusters by the dominant signature (fraction is returned as weight)...")
    if (fit_flag) {
      expo_df <- Signature
    } else {
      expo_df <- get_sig_exposure(Signature, type = "relative")
    }

    sig_names <- colnames(expo_df)[-1]
    common_prefix <- Biobase::lcPrefixC(sig_names)
    mps <- seq_along(sig_names)
    names(mps) <- sig_names

    data <- expo_df %>%
      tidyr::gather(key = "Signature", value = "Exposure", dplyr::starts_with(common_prefix)) %>%
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
      dplyr::mutate(group = as.integer(mps[.data$group])) %>%
      dplyr::arrange(.data$group)

    data$group <- as.character(as.integer(factor(data$group)))
    ztable <- table(data$group, data$enrich_sig)
  } else if (method == "k-means") {
    set.seed(seed = 1024)
    if (fit_flag) {
      expo_df <- Signature
    } else {
      expo_df <- get_sig_exposure(Signature, type = "relative")
    }
    contrib <- expo_df %>%
      as.data.frame() %>%
      tibble::column_to_rownames("sample")
    n_cluster <- ifelse(is.null(n_cluster), ncol(contrib), n_cluster)
    send_info("Running k-means with ", n_cluster, " clusters...")
    contrib.km <- tryCatch(
      kmeans(x = contrib, centers = n_cluster),
      error = function(e) {
        stop("A improper cluster number is set!", call. = FALSE)
      }
    )
    sil_width <- cluster::silhouette(contrib.km$cluster, cluster::daisy(contrib))
    send_info("Generating a table of group and signature contribution (stored in 'map_table' attr):")
    ztable <- contrib.km$centers
    print(ztable)
    send_info("Assigning a group to a signature with the maximum fraction...")
    cluster_df <- as.data.frame(apply(t(ztable), 2, function(x) which(x == max(x))))
    colnames(cluster_df)[1] <- "enrich_sig"
    cluster_df$enrich_sig <- colnames(contrib)[cluster_df$enrich_sig]
    data.table::setDT(x = cluster_df, keep.rownames = TRUE)
    colnames(cluster_df)[1] <- "group"
    data <- as.data.frame(contrib.km$cluster)
    colnames(data)[1] <- "group"
    data$silhouette_width <- signif(sil_width[, "sil_width"], 3)
    data.table::setDT(data, keep.rownames = TRUE)
    colnames(data)[1] <- "sample"
    data$group <- as.character(data$group)
    data <- merge(data, cluster_df, by = "group")
    data.table::setcolorder(data, neworder = c("sample", "group", "silhouette_width", "enrich_sig"))
  }

  data <- data.table::as.data.table(data)
  if (!match_consensus) {
    data <- data[order(as.integer(data$group))]
  }
  send_info("Summarizing...")
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

  send_warning(
    "The 'enrich_sig' column is set to dominant signature in one group, ",
    "please check and make it consistent with biological meaning (correct it by hand if necessary)."
  )
  return(data)
}

utils::globalVariables(
  c("silhouette")
)
