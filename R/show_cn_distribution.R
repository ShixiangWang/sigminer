#' Show Copy Number Distribution either by Length or Chromosome
#'
#' Visually summarize copy number distribution either by copy number segment length
#' or chromosome. Input is a [CopyNumber] object, `genome_build` option will
#' read from `genome_build` slot of object.
#'
#' @param data a [CopyNumber] object.
#' @param rm_normal logical. Whether remove normal copy (i.e. "segVal" equals 2), default is `TRUE`.
#' @param mode either "ld" for distribution by CN length or "cd" for distribution by chromosome.
#' @param fill when `mode` is "cd" and `fill` is `TRUE`, plot percentage instead of count.
#' @param scale_chr logical. If `TRUE`, normalize count to per Megabase unit.
#' @param base_size overall font size.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a `ggplot` object
#' @export
#' @examples
#' # Load copy number object
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Plot distribution
#' p1 <- show_cn_distribution(cn)
#' p1
#' p2 <- show_cn_distribution(cn, mode = "cd")
#' p2
#' p3 <- show_cn_distribution(cn, mode = "cd", fill = TRUE)
#' p3
#' @testexamples
#' expect_s3_class(p1, "ggplot")
#' expect_s3_class(p2, "ggplot")
#' expect_s3_class(p3, "ggplot")
show_cn_distribution <- function(data,
                                 rm_normal = TRUE,
                                 mode = c("ld", "cd"),
                                 fill = FALSE,
                                 scale_chr = TRUE,
                                 base_size = 14) {
  stopifnot(
    is.logical(rm_normal),
    inherits(data, "CopyNumber"),
    is.logical(fill)
  )
  mode <- match.arg(mode)
  ggplot2::theme_set(cowplot::theme_cowplot(font_size = base_size))

  genome_build <- data@genome_build
  data <- data@annotation
  # if remove normal copy number segments
  if (rm_normal) {
    if (!"segVal" %in% colnames(data)) {
      stop("'segVal' must be provided as a column.")
    }
    data <- data[segVal != 2]
  }

  if (mode == "ld") {
    # plot length distribution
    if (!"fraction" %in% colnames(data)) {
      stop("'fraction' must be provided as a column.")
    }
    ggplot(data, aes(x = fraction, y = ..density..)) +
      geom_histogram(bins = 100) +
      labs(
        x = "Length of SCNA\n(normalized to chromosome arms)",
        y = "Percentage\n(as fraction of all SCNAs)"
      )
  } else if (mode == "cd") {
    # plot chr distribution
    if (!all(c("chromosome", "location") %in% colnames(data))) {
      stop("'chromosome', 'location' must be provided as columns.")
    }

    if (is.character(data$chromosome[1])) {
      data$chromosome <- sub("chr", "", data$chromosome, ignore.case = TRUE)
      data <- subset(data, chromosome %in% c(1:22, "X"))
    } else if (is.integer(data$chromosome[1])) {
      data <- subset(data, chromosome %in% c(1:22, "X"))
    }

    # TODO: Find a better way to generate chromosome distribution.

    # only keep chr 1 to 22, X
    data$chromosome <- factor(as.character(data$chromosome),
      levels = c(1:22, "X")
    )

    # only keep p, q, pq
    data$location <- factor(sub("[0-9X]*", "", data$location),
      levels = c("p", "pq", "q")
    )
    # if (sum(!(data$location %in% c("p", "pq", "q"))) > 0){
    #   message("Discarding segments which located in centromere region...")
    #   data = subset(data, location %in% c("p", "pq", "q"))
    # }

    if (scale_chr) {
      chrlen <- get_genome_annotation(
        data_type = "chr_size",
        genome_build = genome_build
      )

      p <- ggplot(data, aes(x = chromosome, fill = location)) +
        geom_bar() +
        xlab("Chromosome")


      q <- ggplot_build(p)$data[[1]][, c("x", "count", "fill")]

      q$x <- factor(q$x, levels = c(1:23), labels = c(1:22, "X"))
      q$fill <- factor(q$fill, levels = c("#F8766D", "#00BA38", "#619CFF"))

      chrlen$chrom <- gsub(
        pattern = "chr",
        replacement = "",
        x = chrlen$chrom
      )
      q <- merge(q, chrlen, by.x = "x", by.y = "chrom")
      q[["count"]] <- 1000000 * (q[["count"]] / q[["size"]])

      if (!fill) {
        ggplot(q, aes(x, y = count, fill = fill)) +
          geom_bar(stat = "identity") +
          scale_fill_discrete(
            name = "location",
            labels = c("p", "pq", "q")
          ) +
          labs(x = "Chromosome", y = "Normalized count (per Mb)")
      } else {
        ggplot(q, aes(x, y = count, fill = fill)) +
          geom_bar(stat = "identity", position = "fill") +
          scale_fill_discrete(
            name = "location",
            labels = c("p", "pq", "q")
          ) +
          labs(x = "Chromosome", y = "Percentage")
      }
    } else {
      # plot
      if (!fill) {
        ggplot(data, aes(x = chromosome, fill = location)) +
          geom_bar() +
          xlab("Chromosome")
      } else {
        ggplot(data, aes(x = chromosome, fill = location)) +
          geom_bar(position = "fill") +
          ylab("Percentage") +
          xlab("Chromosome")
      }
    }
  }
}



utils::globalVariables(
  c(
    "fraction",
    "..density..",
    "location",
    "x",
    "count"
  )
)
