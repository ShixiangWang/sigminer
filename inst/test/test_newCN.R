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

show_cn_features(cn_tally_T, method = "T", nrow = 3)


a = sigminer:::get_cnlist(cn, add_index = T)
getBP10MB_v2(a)
getCNCP_Left(a)
getCNCP_Right(a)
getCNCP_Max(a)

# debug(getOscilation)
getOscilation(a, use_index = T)

test = dplyr::tibble(
  chromosome = c("chr1", "chr1", "chr2", "chr2", "chr2"),
  segVal = c(2, 1, 3, 1, 0)
)

test %>%
dplyr::group_by(.data$chromosome) %>%
  dplyr::mutate(value = abs(c(0L, diff(.data$segVal))))

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


load("~/biodata/DoAbsolute/CN_list.RData")
cn <- read_copynumber(input = CN_list$ACC@data, genome_build = "hg38",
                      seg_cols = c("chromosome", "start", "end", "segVal"))

tally_T <- sig_tally(cn, method = "T")
show_cn_features(tally_T, method = "T", nrow = 3)

cn_list = lapply(CN_list, function(x) {
  read_copynumber(input = x@data, genome_build = "hg38",
                  seg_cols = c("chromosome", "start", "end", "segVal"))
})

tally_list <- lapply(cn_list, function(x) {
  sig_tally(x, method = "T")
})

f_list <- lapply(tally_list, function(x) {
  show_cn_features(x, method = "T", nrow = 3)
})

f_list$ACC
f_list$BLCA
f_list$LAML
f_list$LGG
