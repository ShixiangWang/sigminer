# Source https://github.com/jmonlong/Hippocamplus/blob/master/content/post/2018-06-09-ClusterEqualSize.Rmd

# Todo: support both matrix for data or distance
same_size_clustering <- function(mat, diss = FALSE, clsize = NULL,
                                  algo = c("nnit", "hcbottom", "hclustit", "kmvar"),
                                  method = c(
                                    "maxd", "random", "mind", "elki",
                                    "ward.D", "average", "complete", "single"
                                  )) {
  stopifnot(is.numeric(clsize))

  algo <- match.arg(algo)
  method <- match.arg(method)
  do.call(algo, args = list(mat = mat, diss = diss, clsize = clsize, method = method))
}

nnit <- function(mat, diss = FALSE, clsize = NULL,
                 method = "maxd") {
  clsize.rle <- rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat) / clsize))))
  clsize <- clsize.rle$lengths
  lab <- rep(NA, nrow(mat))
  dmat <- as.matrix(dist(mat))
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

kmvar <- function(mat, diss = FALSE, clsize = NULL,
                  method = "maxd") {
  k <- ceiling(nrow(mat) / clsize)
  km.o <- kmeans(mat, k)
  labs <- rep(NA, nrow(mat))
  centd <- lapply(1:k, function(kk) {
    euc <- t(mat) - km.o$centers[kk, ]
    sqrt(apply(euc, 2, function(x) sum(x^2)))
  })
  centd <- matrix(unlist(centd), ncol = k)
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

hclustit <- function(mat, diss = FALSE, clsize = NULL,
                     method = "ward.D") {
  method <- match.arg(method, choices = c("ward.D", "average", "complete", "single"))
  lab <- rep("l", nrow(mat))
  lab.size <- length(unique(lab))
  while (lab.size > clsize) {
    lab.ii <- which(lab == names(lab.size)[which.max(lab.size)])
    mmat <- mat[lab.ii, ]
    hc.o <- hclust(dist(mmat), method = method)
    lab[lab.ii] <- paste0(lab[lab.ii], "-", cutree(hc.o, 2))
    lab.size <- length(unique(lab))
  }
  lab
}

hcbottom <- function(mat, diss = FALSE, clsize = NULL,
                     method = "ward.D") {
  method <- match.arg(method, choices = c("ward.D", "average", "complete", "single"))
  dmat <- as.matrix(dist(mat))
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
  lab[which(is.na(lab))] <- cpt
  lab
}
