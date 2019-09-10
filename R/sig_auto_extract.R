sig_auto_extract <- function(
                             V0,
                             result_prefix,
                             destdir = tempdir(),
                             K0 = 25,
                             n_run = 10,
                             n_iter = 2e6,
                             tol = 1e-05) {

  # Apply BayesNMF - L1W.L2H for an exponential prior for W and a half-normal prior for H
  for (i in n_run) {
    res <- BayesNMF.L1W.L2H(V0, n_iter, 10, 5, tol, K0, 1)
    save(res, file = file.path(destdir, paste(result_prefix, i, "RData", sep = "."), sep = ""))

    W <- res[[1]] # signature loading
    H <- res[[2]] # activity loading
    index.W <- colSums(W) > 1 # only keep columns in W with non-zero contributions
    W <- W[, index.W]
    H <- H[index.W, ]
    colsum <- colSums(W)
    rowsum <- rowSums(H)

    # By scaling the signature loading matrix has all mutation burdens - each signture (column in W) now represents a number of mutations
    # assigned to each signature.
    for (j in 1:ncol(W)) {
      W[, j] <- W[, j] * rowsum[j]
      H[j, ] <- H[j, ] * colsum[j]
    }
    K <- ncol(W) # number of extracted signatures
    colnames(W) <- paste("W", seq(1:K), sep = "")


    summary.run[i, 1] <- i
    summary.run[i, 2] <- K
    summary.run[i, 3] <- -res[[4]]
  }

  colnames(summary.run) <- c("Run", "K", "posterior")
  summary.run <- data.frame(summary.run)
  # chosen the solution (Run 7) having a maximum posterior at K=12; "L1W.L2H.LEGO96.Biliary.7.RData"
  summary.run <- summary.run[order(summary.run$posterior, decreasing = T), ] ### summary.run is ordered by a posterior
}
