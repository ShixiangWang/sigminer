library(sigminer)
r1 <- sigminer:::bp_extract_signatures(
  t(simulated_catalogs$set1), range = 8:12, cores = 8
)

r1$stats
plot(silhouette ~ signature_number , data = r1$stats, type = "b")
plot(L2_error ~ signature_number , data = r1$stats, type = "b")

sim_results <- purrr::map(
  simulated_catalogs,
  ~ sigminer:::bp_extract_signatures(
    t(.), range = 8:12, cores = 8
  ))

save(sim_results, file = "~/Desktop/sim_results.RData")

sim_results$set1$stats
sim_results$set2$stats
sim_results$set3$stats
sim_results$set4$stats

pheatmap::pheatmap(get_sig_similarity(sim_results$set1$object$K10)$similarity[, c(paste0("COSMIC_", c(1:3, 5:6, 8, 12:13, 17:18)))],
                   cluster_rows = F, cluster_cols = F, display_numbers = T)
pheatmap::pheatmap(get_sig_similarity(sim_results$set2$object$K10)$similarity[, c(paste0("COSMIC_", c(1:3, 5:6, 8, 12:13, 17:18)))],
                   cluster_rows = F, cluster_cols = F, display_numbers = T)
pheatmap::pheatmap(get_sig_similarity(sim_results$set3$object$K10)$similarity[, c(paste0("COSMIC_", c(1:3, 5:6, 8, 12:13, 17:18)))],
                   cluster_rows = F, cluster_cols = F, display_numbers = T)
pheatmap::pheatmap(get_sig_similarity(sim_results$set4$object$K10)$similarity[, c(paste0("COSMIC_", c(1:3, 5:6, 8, 12:13, 17:18)))],
                   cluster_rows = F, cluster_cols = F, display_numbers = T)

# See results without bootstrap
sim_results2 <- purrr::map(
  simulated_catalogs,
  ~ bp_extract_signatures(
    t(.), range = 8:12, cores = 8, n_bootstrap = 0, n_nmf_run = 1000, verbose = TRUE
  ))

sim_results2$set1$stats
sim_results2$set2$stats
sim_results2$set3$stats
sim_results2$set4$stats

