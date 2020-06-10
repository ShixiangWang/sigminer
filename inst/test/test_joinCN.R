load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))
segTabs = segTabs[rev(1:nrow(segTabs))]
cn2 <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE, verbose = TRUE, join_adj_seg = FALSE
)

View(cn@data)

debug(sig_tally)
cn_tally_W <- sig_tally(cn, method = "W")

a <- sigminer:::get_cnlist(cn2)
b <- sigminer:::get_cnlist(segTabs)

identical(a, b)
all.equal(a, b)
str(segTabs$chromosome)

head(a$`TCGA-99-7458-01A-11D-2035-01`)
head(b$`TCGA-99-7458-01A-11D-2035-01`)

b <- sigminer:::get_cnlist(segTabs[, -2])
