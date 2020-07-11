pcawg_cn <- readRDS(file = "../CNSigs/data/pcawg_copynumber.rds")
cn_obj = read_copynumber(pcawg_cn[!Chromosome %in% c("X", "Y")], max_copynumber = 1000L,
                          genome_build = "hg19", complement = FALSE,
                          genome_measure = "called", skip_annotation = TRUE)

extract_seq_dt = function(x) {
  if (inherits(x, "CopyNumber")) {
    x <- x@data
  } else {
    x <- data.table::setDT(x)
  }

  x[, `:=`(
    lenVal = end - start + 1L,
    segVal = ifelse(segVal > 5, 5, segVal) %>% as.integer()  ## Set max value
  )]
  x[, lenVal := cut(lenVal,
                    breaks = c(-Inf, 5e4, 5e5, 5e6, Inf),
                    labels = c("1", "2", "3", "4"),
                    right = FALSE) %>% as.integer()]

  x[, ID := paste(sample, chromosome, sep = "-")]
  x
}


## 4 letters as a group for short, mid, long, long long segments
# A B C D for copy number 0
# E F G H for copy number 1
# I J K L for copy number 2
# M N O P for copy number 3
# Q R S T for copy number 4
# U V W X for copy number 5+
