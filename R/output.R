#' Output Tally Result in Barplots
#'
#' @param x a matrix with row representing components (motifs) and column
#' representing samples.
#' @param result_dir a result directory.
#' @param mut_type one of 'SBS', 'DBS', 'ID' or 'CN'.
#'
#' @return Nothing.
#' @export
output_tally <- function(x, result_dir, mut_type = "SBS") {
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }
  ## output all samples in bar plots
  message("Outputing tally results for ", mut_type)

  if (is.list(x)) {
    x <- x$nmf_matrix %>% t()
  }

  ii <- colSums(x) < 0.01
  if (any(ii)) {
    message(
      "The follow samples dropped due to null catalogue:\n\t",
      paste0(colnames(x)[ii], collapse = ", ")
    )
    x <- x[, !ii, drop = FALSE]
  }

  data.table::fwrite(x %>% as.data.frame() %>% data.table::setDT(keep.rownames = TRUE),
    sep = ",",
    file = file.path(result_dir, paste0(mut_type, "_tally.csv"))
  )

  samps <- colnames(x)

  if (mut_type != "CN") {
    p_total <- show_catalogue(x, style = "cosmic", mode = mut_type, x_label_angle = 90, x_label_vjust = 0.5)
  } else {
    p_total <- show_catalogue(x, style = "cosmic", mode = "copynumber")
  }

  ggsave(file.path(result_dir, paste0(mut_type, "_tally_total.pdf")),
    plot = p_total, width = 12, height = 3,
    limitsize = FALSE
  )

  samp_dir <- file.path(result_dir, paste0(mut_type, "_tally_per_sample"))
  if (!dir.exists(samp_dir)) {
    dir.create(samp_dir, recursive = TRUE)
  }

  for (i in samps) {
    message("Plotting mutation catalog of sample: ", i)
    if (mut_type != "CN") {
      p <- show_catalogue(x,
        style = "cosmic", mode = mut_type,
        samples = i,
        x_label_angle = 90, x_label_vjust = 0.5
      )
    } else {
      p <- show_catalogue(x,
        style = "cosmic", mode = "copynumber",
        normalize = "feature",
        samples = i
      )
    }
    ggsave(file.path(samp_dir, paste0(i, ".pdf")),
      plot = p, width = 12, height = 2, limitsize = FALSE
    )
  }

  return(invisible(NULL))
}

#' Output Signature Results
#'
#' @param sig a `Signature` object.
#' @inheritParams output_tally
#' @inheritParams get_sig_similarity
#'
#' @return Nothing.
#' @export
output_sig <- function(sig, result_dir, mut_type = "SBS", sig_db = mut_type) {
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }
  message("Outputing signature results for ", mut_type, "\n==")

  stopifnot(inherits(sig, "Signature"))

  message("Outputing signature object, signature and exposure matrix.")
  ## Output Signature object
  saveRDS(sig, file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_signatureObj.rds")))

  ## Output data
  if (mut_type != "CN") {
    data.table::fwrite(sig_signature(sig) %>% data.table::as.data.table(keep.rownames = "component"),
      file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_signature.csv"))
    )
  } else {
    data.table::fwrite(sig_signature(sig, normalize = "feature") %>% data.table::as.data.table(keep.rownames = "component"),
      file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_signature.csv"))
    )
  }

  expo_rel <- get_sig_exposure(sig, type = "relative")
  data.table::fwrite(expo_rel,
    file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_relative_exposure.csv"))
  )

  expo_rel_long <- expo_rel %>%
    tidyr::pivot_longer(
      cols = colnames(expo_rel)[-1],
      names_to = "sig",
      values_to = "expo"
    )

  p <- show_group_distribution(expo_rel_long, gvar = "sig", dvar = "expo")
  ggsave(file.path(
    result_dir,
    paste0(
      mut_type, "_",
      attr(sig, "call_method"), "_exposure_distribution_relative.pdf"
    )
  ),
  plot = p, height = 4, width = 1 * sig$K,
  limitsize = FALSE
  )

  expo_abs <- get_sig_exposure(sig)
  data.table::fwrite(expo_abs,
    file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_absolute_exposure.csv"))
  )


  expo_abs_long <- expo_abs %>%
    tidyr::pivot_longer(
      cols = colnames(expo_abs)[-1],
      names_to = "sig",
      values_to = "expo"
    )
  p <- show_group_distribution(expo_abs_long, gvar = "sig", dvar = "expo")
  ggsave(file.path(
    result_dir,
    paste0(
      mut_type, "_",
      attr(sig, "call_method"), "_exposure_distribution_absolute.pdf"
    )
  ),
  plot = p, height = 4, width = 1 * sig$K,
  limitsize = FALSE
  )

  message("Outputing sample clusters based on relative signature exposures.")
  if (sig$K > 1) {
    message("=> Running k-means clustering.")
    grp <- tryCatch(
      get_groups(sig, method = "k-means"),
      error = function(e) {
        message("==> Error detected when running k-means, switch to directly assign samples based on exposure")
        get_groups(sig, method = "exposure")
      }
    )
    data.table::fwrite(grp,
      file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_kmeans_cluster.csv"))
    )
    pdf(
      file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_kmeans_cluster.pdf")),
      width = 1.5 * length(unique(grp$group)), height = 4
    )
    show_groups(grp)
    dev.off()
  } else {
    message("Skip clustering only one signature.")
  }

  ## Output plots
  message("Outputing signature profile plot.")
  if (mut_type != "CN") {
    p <- show_sig_profile(sig, mode = mut_type, style = "cosmic", x_label_angle = 90, x_label_vjust = 0.5)
  } else {
    p <- show_sig_profile(sig, mode = "copynumber", normalize = "feature", style = "cosmic")
  }
  ggsave(file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_signature_profile.pdf")),
    plot = p, width = 12, height = 2 * sig$K,
    limitsize = FALSE
  )

  message("Outputing signature exposure plot.")

  if (ncol(sig$Exposure) < 50) {
    args_expo <- list(sig, hide_samps = FALSE)
  } else {
    args_expo <- list(sig, hide_samps = TRUE, rm_space = TRUE)
  }

  slient_expo <- FALSE
  if (sig$K <= 9) {
    args_expo <- c(args_expo, style = "cosmic")
  } else if (sig$K <= 19) {
    args_expo <- c(args_expo, style = "default")
  } else if (sig$K <= 27) {
    args_expo <- c(args_expo, palette = list(letter_colors))
  } else {
    message("Note: Too many signatures (N>27). Plotting is skipped.")
    message("Please plot it with function `show_sig_exposure()` by your own.")
    slient_expo <- TRUE
  }

  if (!slient_expo) {
    p <- do.call("show_sig_exposure", args_expo)

    cowplot::ggsave2(file.path(
      result_dir, paste0(
        mut_type, "_", attr(sig, "call_method"),
        "_exposure_profile.pdf"
      )
    ),
    plot = p,
    width = 12, height = 6 * (sig$K / 10),
    limitsize = FALSE
    )
  }

  ## Similar analysis and output
  if (!is.null(sig_db)) {
    if (mut_type != "CN") {
      message("Outputing signature similarity analysis results.")
      sim <- get_sig_similarity(sig, sig_db = sig_db)
      data.table::fwrite(sim$similarity %>% data.table::as.data.table(keep.rownames = "sig"),
                         file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_similarity.csv"))
      )
      pheatmap::pheatmap(sim$similarity,
                         cluster_cols = TRUE, cluster_rows = FALSE,
                         filename = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_similarity.pdf")),
                         cellheight = 15, fontsize = 7
      )
      data.table::fwrite(sim$best_match %>% data.table::as.data.table(),
                         file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_COSMIC_best_match.csv"))
      )
    }

    if (mut_type == "SBS" & nchar(rownames(sig$Signature)[1]) == 7) {
      ## Append COSMIC v2 results
      sim <- get_sig_similarity(sig, sig_db = "legacy")
      data.table::fwrite(sim$similarity %>% data.table::as.data.table(keep.rownames = "sig"),
                         file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_legacy_similarity.csv"))
      )
      pheatmap::pheatmap(sim$similarity,
                         cluster_cols = TRUE, cluster_rows = FALSE,
                         filename = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_legacy_similarity.pdf")),
                         cellheight = 15, fontsize = 7
      )
      data.table::fwrite(sim$best_match %>% data.table::as.data.table(),
                         file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_legacy_COSMIC_best_match.csv"))
      )
    }
  } else {
    message("Similarity analysis is skipped.")
  }

  if (attr(sig, "call_method") == "BayesianNMF") {
    data.table::fwrite(
      sig$Raw$summary_run,
      file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_signature_number_survey.csv"))
    )
  }

  if (attr(sig, "call_method") == "SigProfiler") {
    data.table::fwrite(
      sig$Stats$samples,
      file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_sample_stats.csv"))
    )
    data.table::fwrite(
      sig$Stats$signatures,
      file = file.path(result_dir, paste0(mut_type, "_", attr(sig, "call_method"), "_signature_stats.csv"))
    )
  }

  return(invisible(NULL))
}

#' Output Signature Fitting Results
#'
#' @param x result from [sig_fit].
#' @inheritParams output_tally
#' @inheritParams get_sig_similarity
#'
#' @return Nothing.
#' @export
output_fit <- function(x, result_dir, mut_type = "SBS", sig_db = mut_type) {
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }
  message("Outputing signature fitting results for ", mut_type)

  expo <- x$expo
  error <- x$errors %>% data.table::as.data.table(keep.rownames = TRUE)
  colnames(error) <- c("sample", "error")

  rel_expo <- expo %>%
    dplyr::mutate(s = rowSums(.[, -1])) %>%
    dplyr::mutate_at(dplyr::vars(-c("sample", "s")), ~ . / .data$s) %>%
    dplyr::select(-"s") %>%
    data.table::as.data.table()

  data.table::fwrite(expo, file = file.path(result_dir, paste0(mut_type, "_fitting_absolute_exposure.csv")))
  data.table::fwrite(rel_expo, file = file.path(result_dir, paste0(mut_type, "_fitting_relative_exposure.csv")))
  data.table::fwrite(error, file = file.path(result_dir, paste0(mut_type, "_fitting_reconstruction_errors.csv")))

  if (mut_type != "SBS") {
    p1 <- show_sig_fit(expo, palette = NULL, plot_fun = "boxplot") + ggpubr::rotate_x_text()
    p2 <- show_sig_fit(expo, palette = NULL, plot_fun = "violin") + ggpubr::rotate_x_text()

    p3 <- show_sig_fit(rel_expo, palette = NULL, plot_fun = "boxplot") + ggpubr::rotate_x_text()
    p4 <- show_sig_fit(rel_expo, palette = NULL, plot_fun = "violin") + ggpubr::rotate_x_text()
  } else {
    z <- get_sig_db(sig_db)
    sigs <- rownames(z$aetiology)[!grepl("artefact", z$aetiology$aetiology)]
    message("Removed 'Possible sequencing artefact' signatures in plots.")

    p1 <- show_sig_fit(expo, palette = NULL, plot_fun = "boxplot", signatures = sigs) + ggpubr::rotate_x_text()
    p2 <- show_sig_fit(expo, palette = NULL, plot_fun = "violin", signatures = sigs) + ggpubr::rotate_x_text()

    p3 <- show_sig_fit(rel_expo, palette = NULL, plot_fun = "boxplot", signatures = sigs) + ggpubr::rotate_x_text()
    p4 <- show_sig_fit(rel_expo, palette = NULL, plot_fun = "violin", signatures = sigs) + ggpubr::rotate_x_text()
  }

  if (mut_type == "legacy") {
    ## 30 signatures
    width <- 10
    height <- 4
  } else if (mut_type == "SBS") {
    ## 49 signatures
    width <- 15
    height <- 4
  } else if (mut_type == "DBS") {
    ## 11 signatures
    width <- 5
    height <- 4
  } else if (mut_type == "ID") {
    ## 17 signatures
    width <- 7
    height <- 4
  }

  ggsave(file.path(result_dir, paste0(mut_type, "_fitting_absolute_exposure_boxplot.pdf")),
    plot = p1, width = width, height = height,
    limitsize = FALSE
  )
  ggsave(file.path(result_dir, paste0(mut_type, "_fitting_absolute_exposure_violin.pdf")),
    plot = p2, width = width, height = height,
    limitsize = FALSE
  )
  ggsave(file.path(result_dir, paste0(mut_type, "_fitting_relative_exposure_boxplot.pdf")),
    plot = p3, width = width, height = height,
    limitsize = FALSE
  )
  ggsave(file.path(result_dir, paste0(mut_type, "_fitting_relative_exposure_violin.pdf")),
    plot = p4, width = width, height = height,
    limitsize = FALSE
  )
}

#' Output Signature Bootstrap Fitting Results
#'
#' @param x result from [sig_fit_bootstrap_batch].
#' @inheritParams output_tally
#' @inheritParams get_sig_similarity
#' @return Nothing.
#' @export
output_bootstrap <- function(x, result_dir, mut_type = "SBS", sig_db = mut_type) {
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }
  message("Outputing signature bootstrap fitting results for ", mut_type)

  data.table::fwrite(x$expo, file = file.path(result_dir, paste0(mut_type, "_bootstrap_absolute_exposure.csv")))
  data.table::fwrite(x$error, file = file.path(result_dir, paste0(mut_type, "_bootstrap_reconstruction_errors.csv")))
  data.table::fwrite(x$p_val, file = file.path(result_dir, paste0(mut_type, "_bootstrap_p_values_under_different_exposure_cutoffs.csv")))

  if (mut_type != "SBS") {
    p1 <- show_sig_bootstrap_stability(x) + theme(legend.position = "none") + ggpubr::rotate_x_text()
    p2 <- show_sig_bootstrap_exposure(x) + theme(legend.position = "none") + ggpubr::rotate_x_text()
  } else {
    z <- get_sig_db(sig_db)
    sigs <- rownames(z$aetiology)[!grepl("artefact", z$aetiology$aetiology)]
    message("Removed 'Possible sequencing artefact' signatures in plots.")

    p1 <- show_sig_bootstrap_stability(x, signatures = sigs) + theme(legend.position = "none") + ggpubr::rotate_x_text()
    p2 <- show_sig_bootstrap_exposure(x, signatures = sigs) + theme(legend.position = "none") + ggpubr::rotate_x_text()
  }

  if (mut_type == "legacy") {
    ## 30 signatures
    width <- 10
    height <- 4
  } else if (mut_type == "SBS") {
    ## 49 signatures
    width <- 15
    height <- 4
  } else if (mut_type == "DBS") {
    ## 11 signatures
    width <- 5
    height <- 4
  } else if (mut_type == "ID") {
    ## 17 signatures
    width <- 7
    height <- 4
  }

  ggsave(file.path(result_dir, paste0(mut_type, "_bootstrap_signature_instability_boxplot.pdf")),
    plot = p1, width = width, height = height,
    limitsize = FALSE
  )
  ggsave(file.path(result_dir, paste0(mut_type, "_bootstrap_absolute_exposure_boxplot.pdf")),
    plot = p2, width = width, height = height,
    limitsize = FALSE
  )

  samps <- unique(x$expo$sample)
  samp_dir <- file.path(result_dir, paste0(mut_type, "_bootstrap_absolute_exposure_per_sample_boxplot"))
  if (!dir.exists(samp_dir)) {
    dir.create(samp_dir, recursive = TRUE)
  }
  for (i in samps) {
    message("Plotting sample: ", i)
    if (mut_type == "SBS") {
      p <- show_sig_bootstrap_exposure(x, signatures = sigs, sample = i) + theme(legend.position = "none") + ggpubr::rotate_x_text()
    } else {
      p <- show_sig_bootstrap_exposure(x, sample = i) + theme(legend.position = "none") + ggpubr::rotate_x_text()
    }
    ggsave(file.path(samp_dir, paste0(i, ".pdf")),
      plot = p, width = width, height = height,
      limitsize = FALSE
    )
  }
}
