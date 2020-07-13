pcawg_cn <- readRDS(file = "../CNSigs/data/pcawg_copynumber.rds")
cn_obj = read_copynumber(pcawg_cn[!Chromosome %in% c("X", "Y")], max_copynumber = 1000L,
                          genome_build = "hg19", complement = FALSE,
                          genome_measure = "called", skip_annotation = TRUE)

extract_seq_dt = function(x) {
  if (inherits(x, "CopyNumber")) {
    x <- x@data
  } else {
    x <- data.table::setDT(x)
  }

  x[, `:=`(
    lenVal = end - start + 1L,
    segVal = ifelse(segVal > 5, 5, segVal) %>% as.integer()  ## Set max value
  )]
  x[, lenVal := cut(lenVal,
                    breaks = c(-Inf, 5e4, 5e5, 5e6, Inf),
                    labels = c("1", "2", "3", "4"),
                    right = FALSE) %>% as.integer()]

  x[, ID := paste(sample, chromosome, sep = "-")]
  x
}

x <- extract_seq_dt(cn_obj)

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

score_pairwise_strings <- function(x, y, sub_mat) {
  ## For length-1 string "ABCD"
  ## or vector c("A", "B", "C", "D")
  if (length(x) == 1) {
    x <- strsplit(x, "")[[1]]
  }
  if (length(y) == 1) {
    y <- strsplit(y, "")[[1]]
  }
  sub_mat[x, y] %>% diag() %>% sum()
}

get_score_matrix <- function(x, sub_mat, verbose = TRUE) {
  n <- length(x)
  #mat <- matrix(nrow = n, ncol = n, byrow = TRUE)
  mat <- ff::ff(NA_integer_,
                dim = c(n, n), vmode = "byte")  ## Byte from -128 ~ 127
  # Matrix column is faster than row

  i <- j <- 1
  for (i in seq_len(n)) {
    if (verbose) message("Handling sequence: ", x[i])
    j_vals <- vector(mode = "integer", length = i) #ff::ff(NA_integer_, length = i, vmode = "byte")
    for (j in seq_len(i)) {
       j_vals[j] <- score_pairwise_strings(x[i], x[j], sub_mat = sub_mat)
    }
    mat[1:length(j_vals), i] <- j_vals
  }

  rownames(mat) <- colnames(mat) <- x
  return(mat)
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
           by = "ID"]

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

sub_list <- build_sub_matrix()
x$Seqs <- sub_list$map[paste0(x$lenVal, x$segVal)]
x$Seqs

seqs <- extract_sequences(x)

length(seqs$keep)
system.time(
  score_mat <- get_score_matrix(seqs$keep[1:1000], sub_list$mat)
)

system.time(
  score_mat2 <- get_score_matrix(seqs$keep, sub_list$mat, verbose = FALSE)
)

for (i in LETTERS[1:24]) {
  message("Number of seqs starts with ", i, ": #", sum(startsWith(seqs$keep, i)))
}
