extract_seq_dt <- function(x) {
  if (inherits(x, "CopyNumber")) {
    x <- x@data
  } else {
    x <- data.table::setDT(x)
  }

  x[, `:=`(
    lenVal = end - start + 1L,
    segVal = ifelse(segVal > 5, 5, segVal) %>% as.integer() ## Set max value
  )]
  x[, lenVal := cut(lenVal,
    breaks = c(-Inf, 5e4, 5e5, 5e6, Inf),
    labels = c("1", "2", "3", "4"),
    right = FALSE
  ) %>% as.integer()]

  x[, ID := paste(sample, chromosome, sep = "-")]
  x
}

## 4 letters as a group for short, mid, long, long long segments
# A B C D for copy number 0
# E F G H for copy number 1
# I J K L for copy number 2
# M N O P for copy number 3
# Q R S T for copy number 4
# U V W X for copy number 5+

build_sub_matrix <- function() {
  l <- 1:4
  v <- 0:5
  k <- LETTERS[1:24]
  map <- k
  names(map) <- vector_to_combination(1:4, 0:5)

  max_l <- length(l)
  max_v <- length(v)
  pair_mat <- expand.grid(l, v, KEEP.OUT.ATTRS = FALSE) %>% as.matrix()

  score_mat <- pairScoreMatrix(pair_mat, pair_mat, max_l, max_v)
  rownames(score_mat) <- colnames(score_mat) <- k

  return(list(
    map = map,
    mat = score_mat
  ))
}

collapse_shift_seqs <- function(x, len = 5L, step = 1L) {
  if (length(x) <= len) {
    return(paste(x, collapse = ""))
  } else {
    y <- c()
    i <- 0L
    while (i <= (length(x) - len)) {
      y <- c(y, paste(x[(1 + i):(len + i)], collapse = ""))
      i <- i + step
    }
    return(y)
  }
}

extract_sequences <- function(dt, len = 5L, step = 2L, return_dt = FALSE) {
  dt <- dt[, c("ID", "Seqs")]
  dt <- dt[, list(Seqs = collapse_shift_seqs(Seqs, len = len, step = step)),
    by = "ID"
  ]

  if (return_dt) {
    return(dt)
  }

  all_seqs <- unique(dt$Seqs)
  keep <- nchar(all_seqs) >= len

  list(
    keep = sort(all_seqs[keep]),
    drop = sort(all_seqs[!keep])
  )
}

score_pairwise_strings <- function(x, y, sub_mat) {
  ## For length-1 string "ABCD"
  ## or vector c("A", "B", "C", "D")
  if (length(x) == 1) {
    x <- strsplit(x, "")[[1]]
  }
  if (length(y) == 1) {
    y <- strsplit(y, "")[[1]]
  }
  sub_mat[x, y] %>%
    diag() %>%
    sum()
}

get_score_matrix <- function(x, sub_mat, method = c("base", "ff", "bigmemory"), verbose = TRUE) {
  method <- match.arg(method)
  n <- length(x)

  if (method == "base") {
    mat <- matrix(NA_integer_, nrow = n, ncol = n)
  } else if (method == "ff") {
    mat <- ff::ff(NA_integer_,
      dim = c(n, n), vmode = "byte"
    ) ## Byte from -128 ~ 127
  } else {
    options(bigmemory.allow.dimnames = TRUE, bigmemory.typecast.warning = FALSE)
    mat <- bigmemory::big.matrix(n, n, type = "integer")
  }
  # Matrix column is faster than row

  i <- j <- 1
  for (i in seq_len(n)) {
    if (verbose) message("Handling sequence: ", x[i])
    j_vals <- vector(mode = "integer", length = i)
    for (j in seq_len(i)) {
      j_vals[j] <- score_pairwise_strings(x[i], x[j], sub_mat = sub_mat) %>% as.integer()
    }
    mat[seq_along(j_vals), i] <- j_vals
  }

  mat <- mat[]
  ## NOTE the t() operation
  ## cannot just assign upper to lower triangle matrix
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  rownames(mat) <- colnames(mat) <- x
  return(mat)
}

get_score_matrix2 <- function(x, sub_mat, block_size = NULL, verbose = TRUE) {
  if (anyNA(sub_mat)) {
    stop("Input substitution matrix cannot contain 'NA' values!")
  }

  map <- seq_len(24L)
  names(map) <- LETTERS[map]
  map <- map - 1L # to 0 based index

  ## Checking input
  if (any(grepl("[^A-X]", x, ignore.case = FALSE))) {
    stop("The input sequences should contain only A->X, any other letters are invalid.")
  }

  m <- matrix(NA_integer_, ncol = length(x), nrow = nchar(x[1]))

  for (i in seq_len(ncol(m))) {
    s <- unlist(strsplit(x[i], split = ""))
    m[, i] <- map[s] %>% as.integer()
  }
  m <- t(m)

  if (!is.null(block_size)) {
    stopifnot(block_size > 1)
  } else {
    block_size = 1
  }

  y <- getScoreMatrix(m, sub_mat, block_size, verbose)

  if (block_size == 1) {
    colnames(y) <- rownames(y) <- x
  } else {
    colnames(y) <- rownames(y) <- paste0("block", seq_len(nrow(y)))
  }

  return(y)
}

show_seq_shape <- function(x, map = NULL, x_lab = "Estimated segment length", y_lab = "Copy number") {
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required, please install it firstly!")
  }

  if (is.null(map)) {
    map <- LETTERS[1:24]
    names(map) <- vector_to_combination(1:4, 0:5)
  }

  map_df <- data.frame(
    lenVal = strsplit(names(map), split = "") %>% purrr::map_int(~ as.integer(.[1])),
    segVal = strsplit(names(map), split = "") %>% purrr::map_int(~ as.integer(.[2])),
    stringsAsFactors = FALSE
  )
  rownames(map_df) <- as.character(map)

  ## test data
  df <- map_df[unlist(strsplit(x, split = "")), ]
  df$x_end <- cumsum(df$lenVal)
  df$x <- dplyr::lag(df$x_end, default = 0)
  df$color <- ifelse(df$segVal > 2, "red",
    ifelse(df$segVal < 2, "blue",
      "black"
    )
  )

  ggplot(df, aes_string(x = "x", y = "segVal", xend = "x_end", yend = "segVal")) +
    geom_segment(color = df$color) +
    scale_y_continuous(breaks = 0:5, labels = c(0:4, "5+"), limits = c(0, 5)) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(x = x_lab, y = y_lab) +
    cowplot::theme_cowplot()
}

show_seq_logo <- function(x, method = c("prob", "bits"), ncol = NULL, nrow = NULL, ...) {
  method <- match.arg(method)

  if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
    stop("Package 'ggseqlogo' is required, please install it firstly!")
  }

  ## copy from utils.R
  reds <- sapply(list(c(252, 138, 106), c(241, 68, 50), c(188, 25, 26)),
                 FUN = function(x) rgb2hex(x[1], x[2], x[3])) %>% as.character()
  blues <- sapply(list(c(74, 152, 201), c(23, 100, 171)),
                 FUN = function(x) rgb2hex(x[1], x[2], x[3])) %>% as.character()

  cs = ggseqlogo::make_col_scheme(chars = LETTERS[1:24],
                                  groups = c(rep("2 copy DEL", 4),
                                             rep("1 copy DEL", 4),
                                             rep("Normal", 4),
                                             rep("1 copy AMP", 4),
                                             rep("2 copy AMP", 4),
                                             rep("3+ copy AMP", 4)),
                                  cols = c(rep("blue", 4),
                                           rep(blues[1], 4),
                                           rep("black", 4),
                                           rep(reds[1], 4),
                                           rep(reds[2], 4),
                                           rep(reds[3], 4)),
                                  name = "Segment type")

  ggseqlogo::ggseqlogo(x,
                       ncol = ncol,
                       nrow = nrow,
                       method = method,
                       namespace = LETTERS[1:24],
                       col_scheme = cs,
                       ...)
}
