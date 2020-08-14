# Load copy number object
cn <- readRDS("D:/data/pcawg_cn_obj.rds")

cn_tally <- sig_tally(cn, method = "T")

cn_tally$nmf_matrix %>% colnames()
show_catalogue(cn_tally, mode = "copynumber", method = "T", style = "cosmic")

abs_profiles <- sigminer::get_cnlist(cn, add_index = TRUE)
z = cn@data[, grepl(",2,1,2,", paste0(as.integer(segVal), collapse = ",")), by = .(sample, chromosome)]
any(z$V1)

z2 = purrr::map_lgl(abs_profiles, function(x) {
  x2 = x[, grepl(",2,1,2,", paste0(as.integer(segVal),collapse = ",")), by = .(chromosome)]
  any(x2$V1)
})
z2[z2]


cc = paste0(abs_profiles$DO1001$segVal, collapse = ",")
cc

View(abs_profiles$DO1001)

test = list(D1 = abs_profiles$DO1001)

CN_components <- cn_tally$component

