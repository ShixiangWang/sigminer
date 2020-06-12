load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))
segTabs = segTabs[rev(1:nrow(segTabs))]
# segTabs$cc = 1
# segTabs$b = "a"
cn <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE, verbose = TRUE, join_adj_seg = TRUE, use_all = FALSE
)

cn@data

View(cn@data)

cn_tally_T <- sig_tally(cn, method = "T")
cn_tally_T$CN
cn_tally_T$SS

# a <- sigminer:::get_cnlist(cn)
# b <- sigminer:::get_cnlist(segTabs)
#
# identical(a, b)
# all.equal(a, b)
# str(segTabs$chromosome)
#
# head(a$`TCGA-99-7458-01A-11D-2035-01`)
# head(b$`TCGA-99-7458-01A-11D-2035-01`)
#
# b <- sigminer:::get_cnlist(segTabs[, -2])
