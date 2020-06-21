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


load("~/biodata/DoAbsolute/CN_list.RData")
cn <- read_copynumber(input = CN_list$ACC@data, genome_build = "hg38",
                      seg_cols = c("chromosome", "start", "end", "segVal"))

tally_T <- sig_tally(cn, method = "T")
colSums(tally_T$nmf_matrix)

cn_list = lapply(CN_list, function(x) {
  read_copynumber(input = x@data, genome_build = "hg38",
                  seg_cols = c("chromosome", "start", "end", "segVal"))
})

tally_list <- lapply(cn_list, function(x) {
  sig_tally(x, method = "T")
})

sapply(tally_list, function(x) {
  colSums(x$nmf_matrix)
}) %>% t()
