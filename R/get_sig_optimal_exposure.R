# ###############################
# ###### This function is used for the first step of signature attribution to determine an optimal set of signatures out of input signatures (W), given V.
# ###### DEFALUT PARAMETERS: a0=10, tol=1.e-07, phi = 1.50
# ###### INPUTS: V - mutation counts matrix, W - columnwise-normalized signature matrix, Z - binary matrix to indicate a signature availability in each sample.
# ###### OUTPUTS: H - activity-loading matrix
# ###############################
# get_sig_optimal_exposure <- function(nmf_matrix, Signature,
#                                      Z = 1,
#                                      n.iter = 2e5,
#                                      tol = 1e-07,
#                                      phi = 1.5) {
#   V = t(nmf_matrix)
#   W <- Signature$Signature.norm
#
#   a0 <- 10
#   # phi = ? is confusing, follow the orignal code for now
#   K0 <- ncol(W) # of signatures
#
#   # First step: determine a set of optimal signatures best explaining the observed mutations
#   # The automatic relevance determination technique was applied to the H matrix only, while keeping signatures (W) frozen.
#   lambda <- rep(1 + (sqrt((a0 - 1) * (a0 - 2) * mean(V, na.rm = T) / K0)) / (nrow(V) + ncol(V) + a0 + 1) * 2, K0)
#   res0 <- BayesNMF.L1.KL.fixed_W.Z(as.matrix(V), as.matrix(W), as.matrix(H), as.matrix(Z), lambda, n.iter, a0, tol, 1 / phi)
#   H2 <- res0[[2]]
#   colnames(H2) <- colnames(V1)
#   rownames(H2) <- colnames(W0)
#
#   # Second step: determine a sample-level signature attribution using selected sigantures in the first step.
#   index.H2 <- rowSums(H2) > 1 ### identify only active signatures
#   Z2 <- Z1
#   Z2[!index.H2, ] <- 0 ### only selected signatures in the first step are allowed + the original contraints on the signature availability from Z1.
#   for (j in 1:ncol(H2)) {
#     tmpH <- rep(0, ncol(W0))
#     if (sum(V1[, j]) >= 5) {
#       lambda <- 1 + (sqrt((a0 - 1) * (a0 - 2) * mean(V1, na.rm = T) / K0)) / (nrow(V1) + ncol(V1) + a0 + 1) * 2
#       res <- BayesNMF.L1.KL.fixed_W.Z.sample(
#         as.matrix(V1[, j]), W,
#         as.matrix(H2[, j]), as.matrix(Z2[, j]),
#         lambda, n.iter, a0, tol, 1
#       )
#       tmpH <- res[[2]]
#     }
#     if (j == 1) {
#       H3 <- tmpH
#     } else {
#       H3 <- cbind(H3, tmpH)
#       cat(j, "\n")
#     }
#   }
#   colnames(H3) <- colnames(V1)
#   rownames(H3) <- colnames(W0)
#   H3
# }
