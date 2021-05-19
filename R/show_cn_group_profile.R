#' Show Summary Copy Number Profile for Sample Groups
#'
#' @inheritParams show_cn_freq_circos
#' @inheritParams show_cn_profile
#' @param fill_area default is `TRUE`, fill area with colors.
#' @param force_y_limit default is `TRUE`, force multiple plots
#' @param highlight_genes gene list to highlight.
#' have same y ranges. You can also set a length-2 numeric value.
#' @param repel if `TRUE` (default is `FALSE`), repel highlight genes to
#' avoid overlap.
#'
#' @return a (list of) `ggplot` object.
#' @export
#'
#' @examples
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' p1 <- show_cn_group_profile(cn)
#' p1
#' \donttest{
#' ss <- unique(cn@data$sample)
#' p2 <- show_cn_group_profile(cn, groups = list(a = ss[1:5], b = ss[6:10]))
#' p2
#' p3 <- show_cn_group_profile(cn,
#'   groups = list(g1 = ss[1:5], g2 = ss[6:10]),
#'   force_y_limit = c(-1, 1), nrow = 2
#' )
#' p3
#'
#' ## Set custom cutoff for custom data
#' data <- cn@data
#' data$segVal <- data$segVal - 2L
#' p4 <- show_cn_group_profile(data,
#'   groups = list(g1 = ss[1:5], g2 = ss[6:10]),
#'   force_y_limit = c(-1, 1), nrow = 2,
#'   cutoff = c(0, 0)
#' )
#' p4
#'
#' ## Add highlight gene
#' p5 <- show_cn_group_profile(cn, highlight_genes = c("TP53", "EGFR"))
#' p5
#' }
#' @testexamples
#' expect_s3_class(p1, "ggplot")
#' expect_s3_class(p2, "ggplot")
#' expect_s3_class(p3, "ggplot")
#' expect_s3_class(p4, "ggplot")
#' expect_s3_class(p5, "ggplot")
show_cn_group_profile <- function(data,
                                  groups = NULL,
                                  fill_area = TRUE,
                                  cols = NULL,
                                  chrs = paste0("chr", c(1:22, "X")),
                                  genome_build = c("hg19", "hg38", "mm10", "mm9"),
                                  cutoff = 2L,
                                  resolution_factor = 1L,
                                  force_y_limit = TRUE,
                                  highlight_genes = NULL,
                                  repel = FALSE,
                                  nrow = NULL, ncol = NULL,
                                  return_plotlist = FALSE) {
  stopifnot(is.data.frame(data) | inherits(data, "CopyNumber"))
  if (is.data.frame(data)) {
    nc_cols <- c("chromosome", "start", "end", "segVal", "sample")
    if (!all(nc_cols %in% colnames(data))) {
      stop("Invalid input, it must contain columns: ", paste(nc_cols, collapse = " "))
    }
  }

  if (is.null(cols)) {
    cols <- c("#FF000080", "#0000FF80")
  }

  genome_build <- match.arg(genome_build)
  if (inherits(data, "CopyNumber")) {
    genome_build <- data@genome_build
    data <- data@data
  } else {
    data <- data.table::as.data.table(data)
  }

  # The 'Groups' can be either a list of name index
  # or the column name
  if (is.list(groups)) {
    if (length(purrr::reduce(groups, intersect)) != 0) {
      stop("Duplicated sample names in group list!", call. = FALSE)
    }
    samples <- purrr::reduce(groups, c)
    data <- data[data$sample %in% samples]
    grp_df <- tibble::enframe(groups) %>%
      tidyr::unnest("value") %>%
      purrr::set_names(c("grp_name", "sample")) %>%
      data.table::as.data.table()

    data <- merge(data, grp_df, by = "sample")
  } else if (is.character(groups)) {
    if (length(groups) != 1) {
      stop("When 'groups' is character type, it represents one column name.",
        call. = FALSE
      )
    }
    samples <- unique(data$sample)
    colnames(data)[colnames(data) == groups] <- "grp_name"
  } else {
    samples <- unique(data$sample)
    data$grp_name <- "1"
  }

  data$sample <- factor(data$sample, levels = samples)

  data$chromosome <- ifelse(startsWith(data$chromosome, prefix = "chr"),
    data$chromosome,
    paste0("chr", data$chromosome)
  )

  data <- data[data$chromosome %in% chrs]

  if (length(cutoff) == 1) {
    cutoff <- c(cutoff, cutoff)
  }
  ## Construct CNV frequency data
  ## for AMP/DEL and each group
  data_freq <- data[, get_cn_freq_table(.SD,
    genome_build = genome_build,
    cutoff = cutoff,
    resolution_factor = resolution_factor
  ), by = "grp_name"]
  data.table::setcolorder(data_freq, c(
    "chromosome", "start", "end",
    "freq_AMP", "freq_DEL"
  ))

  if (!is.null(highlight_genes)) {
    gene_file <- switch(genome_build,
      mm10 = file.path(
        system.file("extdata", package = "sigminer"),
        "mouse_mm10_gene_info.rds"
      ),
      mm9 = file.path(
        system.file("extdata", package = "sigminer"),
        "mouse_mm9_gene_info.rds"
      ),
      file.path(
        system.file("extdata", package = "sigminer"),
        paste0("human_", genome_build, "_gene_info.rds")
      )
    )
    ok <- TRUE
    if (!file.exists(gene_file)) ok <- query_remote_data(basename(gene_file))
    if (!ok) {
      return(invisible(NULL))
    }
    gene_dt <- readRDS(gene_file)
    gene_dt <- gene_dt[gene_dt$gene_name %in% highlight_genes][
      , c("chrom", "start", "end", "gene_name")
    ]

    if (nrow(gene_dt) < 1) {
      stop(paste("No gene matched in", genome_build))
    }

    colnames(gene_dt) <- c("chromosome", "i.start", "i.end", "gene_name")
    data_freq <- merge(data_freq, gene_dt, by = "chromosome", all.x = TRUE)
    data_freq$highlight <- data.table::fifelse(
      data_freq$i.end <= data_freq$end & data_freq$i.start >= data_freq$start,
      data_freq$gene_name, NA_character_
    )
    data_freq$i.start <- data_freq$i.end <- data_freq$gene_name <- NULL
  } else {
    data_freq$highlight <- NA_character_
  }

  coord_df <- build_chrom_coordinate(genome_build, chrs) %>%
    dplyr::mutate(labels = sub("chr", "", .data$chrom))

  merge_df <- data_freq %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(coord_df, by = c("chromosome" = "chrom")) %>%
    dplyr::mutate(
      start = .data$x_start + .data$start - 1,
      end = .data$x_start + .data$end
    )
  merge_df$freq_DEL <- -merge_df$freq_DEL

  grp_data <- dplyr::group_split(merge_df, .data$grp_name)

  if (length(grp_data) == 1) {
    gg <- plot_cn_summary(grp_data[[1]], coord_df, fill_area = fill_area, cols = cols, repel = repel)
  } else {
    gglist <- purrr::map(
      grp_data,
      ~ plot_cn_summary(., coord_df, fill_area = fill_area, cols = cols, repel = repel) +
        ggplot2::labs(title = .$grp_name[1])
    )

    if (isTRUE(force_y_limit)) {
      ylim <- range(c(range(merge_df$freq_AMP), range(merge_df$freq_DEL)))
    } else if (is.numeric(force_y_limit)) {
      ylim <- force_y_limit
    }

    gglist <- purrr::map(gglist, ~ . + ggplot2::ylim(ylim))

    if (return_plotlist) {
      return(gglist)
    }
    gg <- cowplot::plot_grid(
      plotlist = gglist, align = "hv",
      nrow = nrow, ncol = ncol
    )
  }

  return(gg)
}


plot_cn_summary <- function(plot_df, coord_df, fill_area = TRUE, cols = c("red", "blue"),
                            repel = FALSE) {
  h_df <- na.omit(plot_df)

  plot_df <- dplyr::bind_rows(
    plot_df %>%
      dplyr::select(c("start", "freq_AMP", "freq_DEL")) %>%
      dplyr::rename(x = .data$start),
    plot_df %>%
      dplyr::select(c("end", "freq_AMP", "freq_DEL")) %>%
      dplyr::rename(x = .data$end) %>%
      dplyr::mutate(x = .data$x - 1)
  ) %>%
    na.omit()

  data_amp <- plot_df %>%
    dplyr::select(c("x", "freq_AMP")) %>%
    purrr::set_names(c("x", "freq"))
  data_del <- plot_df %>%
    dplyr::select(c("x", "freq_DEL")) %>%
    purrr::set_names(c("x", "freq"))

  if (!fill_area) {
    fill_amp_col <- fill_del_col <- "white"
  } else {
    fill_amp_col <- cols[1]
    fill_del_col <- cols[2]
  }

  p <- ggplot() +
    geom_area(aes_string(x = "x", y = "freq"),
      fill = fill_amp_col, color = cols[1], data = data_amp
    ) +
    geom_area(aes_string(x = "x", y = "freq"),
      fill = fill_del_col, color = cols[2], data = data_del
    )

  if (nrow(h_df) > 0) {
    text_fun <- if (repel) ggrepel::geom_text_repel else geom_text
    p <- p + text_fun(
      aes(
        x = (.data$start + .data$end) / 2,
        y = max(data_amp$freq) - 0.02,
        label = .data$highlight
      ),
      size = 2, color = "black",
      data = h_df
    )
  }

  p <- p +
    geom_hline(yintercept = 0) +
    geom_vline(aes(xintercept = .data$x_start),
      alpha = 0.5,
      linetype = "dotted", data = coord_df
    ) +
    geom_vline(
      xintercept = coord_df$x_end[nrow(coord_df)],
      alpha = 0.5,
      linetype = "dotted"
    ) +
    scale_x_continuous(breaks = coord_df$lab_loc, labels = coord_df$labels) +
    labs(x = "Chromosome", y = "Variation frequency") +
    cowplot::theme_cowplot() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 60, hjust = 1, size = 9)
    )
}
