# https://github.com/jmonlong/Hippocamplus/blob/master/content/post/2018-06-09-ClusterEqualSize.Rmd

#' Same Size Clustering
#'
#' This is a wrapper for several implementation that classify samples into
#' same size clusters, the details please see [this blog](http://jmonlong.github.io/Hippocamplus/2018/06/09/cluster-same-size/).
#' The source code is modified based on code from the blog.
#'
#' @param mat a data/distance matrix.
#' @param diss if `TRUE`, treat `mat` as a distance matrix.
#' @param clsize integer, number of sample within a cluster.
#' @param algo algorithm.
#' @param method method.
#'
#' @return a vector.
#' @export
#'
#' @examples
#' set.seed(1234L)
#' x <- rbind(
#'   matrix(rnorm(100, sd = 0.3), ncol = 2),
#'   matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2)
#' )
#' colnames(x) <- c("x", "y")
#'
#' y1 <- same_size_clustering(x, clsize = 10)
#' y11 <- same_size_clustering(as.matrix(dist(x)), clsize = 10, diss = TRUE)
#'
#' y2 <- same_size_clustering(x, clsize = 10, algo = "hcbottom", method = "ward.D")
#'
#' y3 <- same_size_clustering(x, clsize = 10, algo = "kmvar")
#' y33 <- same_size_clustering(as.matrix(dist(x)), clsize = 10, algo = "kmvar", diss = TRUE)
#' @testexamples
#' expect_length(y1, 100L)
#' expect_length(y11, 100L)
#' expect_length(y2, 100L)
#' expect_length(y3, 100L)
#' expect_length(y33, 100L)
same_size_clustering <- function(mat, diss = FALSE, clsize = NULL,
                                 algo = c("nnit", "hcbottom", "kmvar"),
                                 method = c(
                                   "maxd", "random", "mind", "elki",
                                   "ward.D", "average", "complete", "single"
                                 )) {
  stopifnot(is.numeric(clsize))

  algo <- match.arg(algo)
  method <- match.arg(method)
  do.call(algo, args = list(mat = mat, diss = diss, clsize = clsize, method = method))
}

nnit <- function(mat,
                 clsize = NULL,
                 diss = FALSE,
                 method = "maxd") {
  stopifnot(is.logical(diss))

  clsize.rle <- rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat) / clsize))))
  clsize <- clsize.rle$lengths
  lab <- rep(NA, nrow(mat))
  if (isFALSE(diss)) {
    dmat <- as.matrix(dist(mat))
  } else {
    dmat <- mat
  }
  cpt <- 1
  while (sum(is.na(lab)) > 0) {
    lab.ii <- which(is.na(lab))
    dmat.m <- dmat[lab.ii, lab.ii]
    ii <- switch(
      method,
      maxd = which.max(rowSums(dmat.m)),
      mind = which.min(rowSums(dmat.m)),
      random = sample.int(nrow(dmat.m), 1),
      stop("unsupported method in 'nnit'!")
    )
    lab.m <- rep(NA, length(lab.ii))
    lab.m[head(order(dmat.m[ii, ]), clsize[cpt])] <- cpt
    lab[lab.ii] <- lab.m
    cpt <- cpt + 1
  }
  if (any(is.na(lab))) {
    lab[which(is.na(lab))] <- cpt
  }
  lab
}

kmvar <- function(mat,
                  clsize = NULL,
                  diss = FALSE,
                  method = "maxd") {
  stopifnot(is.logical(diss))

  k <- ceiling(nrow(mat) / clsize)
  if (isFALSE(diss)) {
    km.o <- kmeans(mat, k)
    # distance to centers
    centd <- lapply(1:k, function(kk) {
      euc <- t(mat) - km.o$centers[kk, ]
      sqrt(apply(euc, 2, function(x) sum(x^2)))
    })
    centd <- matrix(unlist(centd), ncol = k)
  } else {
    message("PAM algorithm is applied when input distance matrix.")
    pam.o <- cluster::pam(mat, k, diss = TRUE)
    # medoids
    # distance to medoids
    centd <- mat[, pam.o$id.med, drop = FALSE]
  }

  labs <- rep(NA, nrow(mat))
  clsizes <- rep(0, k)

  ptord <- switch(
    method,
    maxd = order(-apply(centd, 1, max)),
    mind = order(apply(centd, 1, min)),
    random = sample.int(nrow(mat)),
    elki = order(apply(centd, 1, min) - apply(centd, 1, max)),
    stop("unsupported method in 'kmvar'!")
  )

  for (ii in ptord) {
    bestcl <- which.max(centd[ii, ])
    labs[ii] <- bestcl
    clsizes[bestcl] <- clsizes[bestcl] + 1
    if (clsizes[bestcl] >= clsize) {
      centd[, bestcl] <- NA
    }
  }
  return(labs)
}

hcbottom <- function(mat,
                     clsize = NULL,
                     diss = FALSE,
                     method = "ward.D") {
  stopifnot(is.logical(diss))

  method <- match.arg(method, choices = c("ward.D", "average", "complete", "single"))
  if (isFALSE(diss)) {
    dmat <- as.matrix(dist(mat))
  } else {
    dmat <- mat
  }
  clsize.rle <- rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat) / clsize))))
  clsizes <- clsize.rle$lengths
  cpt <- 1
  lab <- rep(NA, nrow(mat))
  for (clss in clsizes[-1]) {
    lab.ii <- which(is.na(lab))
    hc.o <- hclust(as.dist(dmat[lab.ii, lab.ii]), method = method)
    clt <- 0
    ct <- length(lab.ii) - clss
    while (max(clt) < clss) {
      cls <- cutree(hc.o, ct)
      clt <- table(cls)
      ct <- ct - 1
    }
    cl.sel <- which(cls == as.numeric(names(clt)[which.max(clt)]))
    lab[lab.ii[head(cl.sel, clss)]] <- cpt
    cpt <- cpt + 1
  }
  lab[is.na(lab)] <- cpt
  lab
}
