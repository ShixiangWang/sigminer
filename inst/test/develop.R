##########################################
### Store temp thoughts and functions ####
##########################################

# set order and names for signature profile
# Load copy number signature
load(system.file("extdata", "toy_copynumber_signature.RData",
                 package = "sigminer", mustWork = TRUE
))
draw_sig_profile(res$nmfObj, sig_names = c("CN.Sig.2", "CN.Sig.1"), sig_orders = c("CN.Sig.2", "CN.Sig.1"))


# feature enrichment and plot (for categorical variables)
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
tt = clinicalEnrichment(laml, 'FAB_classification')


am_data = sig_group_enrichment(test_df)
plotEnrichmentResults(am_data, showTitle = TRUE)


# Load copy number signature
load(system.file("extdata", "toy_copynumber_signature.RData",
                 package = "sigminer", mustWork = TRUE
))
# Assign samples to clusters
subtypes <- sig_assign_samples(res$nmfObj, type = "samples")

set.seed(1234)
# Add custom groups
subtypes$new_group <- sample(c("1", "2", "3", "4", NA), size = nrow(subtypes), replace = TRUE)
# Summarize subtypes
subtypes.sum2 <- sig_summarize_subtypes(subtypes[, -1],
                                       col_subtype = "Subgroup",
                                       cols_to_summary = colnames(subtypes[, -1])[c(-1, -2)],
                                       type = c("co", "ca"), verbose = TRUE, NAs = "NA"
)
tt = draw_subtypes_comparison(subtypes.sum2, label="p.signif")
tt$ca$new_group
tt$co$prob
