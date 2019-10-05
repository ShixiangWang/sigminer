############################################################################################
############################################################################################
#### Copyright (c) 2017, Broad Institute
#### Redistribution and use in source and binary forms, with or without
#### modification, are permitted provided that the following conditions are
#### met:
####     Redistributions of source code must retain the above copyright
####     notice, this list of conditions and the following disclaimer.
####     Redistributions in binary form must reproduce the above copyright
####     notice, this list of conditions and the following disclaimer in
####     the documentation and/or other materials provided with the
####     distribution.
####     Neither the name of the Broad Institute nor the names of its
####     contributors may be used to endorse or promote products derived
####     from this software without specific prior written permission.
#### THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#### "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#### LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#### A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#### HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#### SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#### LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#### DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#### THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#### (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#### OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################################################################################
############################################################################################


# Functions listed in this file are obtained from reference
# Alexandrov, Ludmil, et al. "The repertoire of mutational signatures in human cancer." BioRxiv (2018): 322859.
# or Broad Institute.
# For details on the implementation
# see J Kim, Mouw K, P Polak et al, Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors
# Nat. Genet. (2016) DOI: 10.1038/ng.3557
# For details on the original algorithms
# see Tan, V.Y. & Fe<U+0301>votte, C. Automatic relevance determination in nonnegative matrix factorization with the beta-divergence.
# IEEE Trans. Pattern Anal. Mach. Intell. 35, 1592–1605 (2013).
#
# They are modified by me for better used in sigminer


###############################
###### This function is used in both PRIMARY and SECONDARY steps of the signature extraction in 2780 PCAWG samples.
###### Bayesian non-negative matrix factoriztion algorithm with an exponential prior for W and a half-normal prior for H.
###### DEFALUT PARAMETERS: a0=10, b0=5, phi=1, n.iter=2000000, tol=1.e-05
###### INPUT: V (mutation counts matrix), K0 (maximum num of signatures usually 96)
###### OUTPUT: W (signature-loading matrix), H (activity-loading matrix)
###############################

BayesNMF.L1W.L2H <- function(V0, n.iter, a0, b0, tol, K, K0, phi) {
  eps <- 1.e-50
  del <- 1.0
  V <- V0 - min(V0)
  N <- dim(V)[1]
  M <- dim(V)[2]
  W <- matrix(runif(N * K) * sqrt(mean(V)), ncol = K)
  H <- matrix(runif(M * K) * sqrt(mean(V)), ncol = M)
  V.ap <- W %*% H + eps
  I_NM <- array(1, dim = c(N, M))
  I_NK <- array(1, dim = c(N, K))

  C <- N + M / 2 + a0 - 1
  beta <- C / (colSums(W) + 0.5 * rowSums(H * H) + b0)
  beta.bound <- C / b0
  beta.cut <- beta.bound / 1.25

  n.beta <- list()
  n.beta[[1]] <- beta
  iter <- 2
  while (del >= tol & iter < n.iter) {
    B <- diag(beta)
    H <- H * (t(W) %*% (V / V.ap) / (t(W) %*% I_NM + B %*% H + eps))
    V.ap <- W %*% H + eps
    W <- W * ((V / V.ap) %*% t(H) / (I_NM %*% t(H) + I_NK %*% B + eps))
    beta <- C / (colSums(W) + 0.5 * rowSums(H * H) + b0)
    V.ap <- W %*% H + eps
    n.beta[[iter]] <- beta
    error <- sum((V - V.ap)^2)
    del <- max(abs(beta - n.beta[[iter - 1]]) / n.beta[[iter - 1]])
    like <- sum(V * log((V + eps) / (V.ap + eps)) + V.ap - V)
    evid <- like + sum((colSums(W) + 0.5 * rowSums(H * H) + b0) * beta - C * log(beta))
    if (iter %% 10000 == 0) {
      cat(iter, evid, like, error, del, sum(colSums(W) > 1.e-05), sum(beta <= beta.cut), "
")
    }
    iter <- iter + 1
  }
  lambda <- 1 / beta
  return(list(W, H, like, evid, lambda, error))
}

## Bayesian NMF algorithm with exponential priors for W and H
BayesNMF.L1KL <- function(V0, n.iter, a0, tol, K, K0, phi) {
  eps <- 1.e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[, active_nodes]
  V <- V0 - min(V0) + eps
  Vmax <- max(V)
  N <- dim(V)[1]
  M <- dim(V)[2]
  W <- matrix(runif(N * K) * sqrt(Vmax), ncol = K)
  H <- matrix(runif(M * K) * sqrt(Vmax), ncol = M)
  V.ap <- W %*% H + eps

  C <- N + M + a0 + 1
  b0 <- sqrt((a0 - 1) * (a0 - 2) * mean(V, na.rm = T) / K0)
  lambda.bound <- b0 / C
  lambda <- (colSums(W) + rowSums(H) + b0) / C
  lambda.cut <- 1.5 * lambda.bound

  n.lambda <- list()
  n.lambda[[1]] <- lambda

  iter <- 2
  while ((del >= tol) & (iter < n.iter)) {
    H <- H * (t(W) %*% (V / V.ap)) / (matrix(rep(colSums(W) + phi / lambda, M), ncol = M) + eps)
    V.ap <- W %*% H + eps
    W <- W * ((V / V.ap) %*% t(H)) / t(matrix(rep(rowSums(H) + phi / lambda, N), ncol = N) + eps)
    V.ap <- W %*% H + eps
    lambda <- (colSums(W) + rowSums(H) + b0) / C
    del <- max(abs(lambda - n.lambda[[iter - 1]]) / n.lambda[[iter - 1]])
    like <- sum(V * log(V / V.ap) + V.ap - V)
    evid <- like + phi * sum((colSums(W) + rowSums(H) + b0) / lambda + C * log(lambda))
    n.lambda[[iter]] <- lambda
    error <- sum((V - V.ap)^2)
    if (iter %% 10000 == 0) {
      cat(iter, evid, like, error, del, sum(colSums(W) != 0), sum(lambda >= lambda.cut), "
")
    }
    iter <- iter + 1
  }
  return(list(W, H, like, evid, lambda, error))
}

## Bayesian NMF algorithm with hlaf-normal priors for W and H
BayesNMF.L2KL <- function(V0, n.iter, a0, tol, K, K0, phi) {
  eps <- 1.e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[, active_nodes]
  V <- V0 - min(V0)
  Vmax <- mean(V)
  N <- dim(V)[1]
  M <- dim(V)[2]

  W <- matrix(runif(N * K) * Vmax, ncol = K)
  H <- matrix(runif(M * K) * Vmax, ncol = M)
  V.ap <- W %*% H + eps

  C <- (N + M) / 2 + a0 + 1
  b0 <- 3.14 * (a0 - 1) * mean(V) / (2 * K0)
  lambda <- (0.5 * colSums(W * W) + 0.5 * rowSums(H * H) + b0) / C
  lambda.cut <- b0 / C * 1.25

  I <- array(1, dim = c(N, M))

  n.lambda <- list()
  n.lambda[[1]] <- lambda
  iter <- 2
  while ((del >= tol) & (iter < n.iter)) {
    H <- H * (t(W) %*% (V / V.ap) / (t(W) %*% I + phi * H / matrix(rep(lambda, M), ncol = M) + eps))^0.5
    V.ap <- W %*% H + eps
    W <- W * ((V / V.ap) %*% t(H) / (I %*% t(H) + phi * W / t(matrix(rep(lambda, N), ncol = N)) + eps))^0.5
    lambda <- (0.5 * colSums(W * W) + 0.5 * rowSums(H * H) + b0) / C
    V.ap <- W %*% H + eps
    del <- max(abs(lambda - n.lambda[[iter - 1]]) / n.lambda[[iter - 1]])
    like <- sum(V * log((V + eps) / (V.ap + eps)) + V.ap - V)
    evid <- like + phi * sum((0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / lambda + C * log(lambda))
    n.lambda[[iter]] <- lambda
    error <- sum((V - V.ap)^2)
    if (iter %% 10000 == 0) {
      cat(iter, evid, like, error, del, sum(colSums(W) != 0), sum(lambda >= lambda.cut), "
")
    }
    iter <- iter + 1
  }
  return(list(W, H, like, evid, lambda, error))
}


###############################
###### This function is used for the first step of signature attribution to determine an optimal set of signatures out of input signatures (W), given V.
###### DEFALUT PARAMETERS: a0=10, tol=1.e-07, phi = 1.50
###### INPUTS: V - mutation counts matrix, W - columnwise-normalized signature matrix, Z - binary matrix to indicate a signature availability in each sample.
###### OUTPUTS: H - activity-loading matrix
###############################
BayesNMF.L1.KL.fixed_W.Z <- function(V, W, H, Z, lambda, n.iter, a0, tol, phi) {
  eps <- 1.e-50
  del <- 1.0
  N <- dim(V)[1]
  M <- dim(V)[2]
  V <- V + eps
  K <- ncol(W)
  K0 <- K
  V.ap <- W %*% H + eps
  C <- N + M + a0 + 1
  b0 <- sqrt((a0 - 1) * (a0 - 2) * mean(V, na.rm = T) / K0)
  lambda.bound <- b0 / C
  lambda.cut <- 1.1 * lambda.bound

  n.lambda <- list()
  n.lambda[[1]] <- lambda
  iter <- 2
  while (del >= tol & iter < n.iter) {
    H <- H * (t(W) %*% (V / V.ap)) / (matrix(rep(colSums(W) + phi / lambda, M), ncol = M) + eps)
    H <- H * Z
    V.ap <- W %*% H + eps
    lambda <- (colSums(W) + rowSums(H) + b0) / C
    del <- max(abs(lambda - n.lambda[[iter - 1]]) / n.lambda[[iter - 1]])
    n.lambda[[iter]] <- lambda
    if (iter %% 100 == 0) {
      error <- sum((V - V.ap)^2)
      like <- sum(V * log((V + 1) / (V.ap + 1)) + V.ap - V)
      evid <- like / phi + sum((colSums(W) + rowSums(H) + b0) / lambda + C * log(lambda))
      cat(iter, evid, like, error, del, sum(rowSums(H) != 0), sum(lambda >= lambda.cut), "
")
    }
    iter <- iter + 1
  }
  error <- sum((V - V.ap)^2)
  like <- sum(V * log((V + 1) / (V.ap + 1)) + V.ap - V)
  evid <- like / phi + sum((colSums(W) + rowSums(H) + b0) / lambda + C * log(lambda))
  cat("***********************************************************", "
")
  return(list(W, H, like, evid, lambda, error))
}

###############################
###### This function is used for the second step of signature attribution to determine a sample-level attribution using selected signatures in the first step.
###### DEFALUT PARAMETERS: a0=10, tol=1.e-07, phi=1.0
###### INPUTS: V - mutation counts matrix, W - columnwise-normalized signature matrix, Z - binary matrix to constraint a signature availability.
###### OUTPUTS: H - activity-loading matrix
###############################
BayesNMF.L1.KL.fixed_W.Z.sample <- function(V, W, H, Z, lambda, n.iter, a0, tol, phi) {
  eps <- 1.e-50
  del <- 1.0
  N <- dim(V)[1]
  M <- dim(V)[2]
  V <- V + eps
  K <- ncol(W)
  K0 <- K
  V.ap <- W %*% H + eps
  C <- N + M + a0 + 1
  b0 <- sqrt((a0 - 1) * (a0 - 2) * mean(V, na.rm = T) / K0)
  lambda.bound <- b0 / C
  lambda.cut <- 1.1 * lambda.bound

  n.lambda <- list()
  n.lambda[[1]] <- lambda
  iter <- 2

  Z0 <- Z
  while (del >= tol & iter < n.iter) {
    H <- H * (t(W) %*% (V / V.ap)) / (matrix(rep(colSums(W) + 0 * phi / lambda, M), ncol = M) + eps)
    H <- H * (Z0)
    V.ap <- W %*% H + eps
    lambda <- (colSums(W) + rowSums(H) + b0) / C
    del <- max(abs(lambda - n.lambda[[iter - 1]]) / n.lambda[[iter - 1]])
    n.lambda[[iter]] <- lambda
    if (iter %% 100 == 0) {
      error <- sum((V - V.ap)^2)
      like <- sum(V * log((V + 1) / (V.ap + 1)) + V.ap - V)
      evid <- like / phi + sum((colSums(W) + rowSums(H) + b0) / lambda + C * log(lambda))
      cat(iter, evid, like, error, del, sum(rowSums(H) != 0), sum(lambda >= lambda.cut), "
")
      res <- list(W, H, like, evid, lambda, error)
    }
    iter <- iter + 1
  }
  error <- sum((V - V.ap)^2)
  like <- sum(V * log((V + 1) / (V.ap + 1)) + V.ap - V)
  evid <- like / phi + sum((colSums(W) + rowSums(H) + b0) / lambda + C * log(lambda))
  cat("***********************************************************", "
")
  return(list(W, H, like, evid, lambda, error))
}
