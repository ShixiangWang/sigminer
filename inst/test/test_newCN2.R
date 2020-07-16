pcawg_cn <- readRDS(file = "../CNSigs/data/pcawg_copynumber.rds")
cn_obj = read_copynumber(pcawg_cn[!Chromosome %in% c("X", "Y")], max_copynumber = 1000L,
                          genome_build = "hg19", complement = FALSE,
                          genome_measure = "called", skip_annotation = TRUE)

x <- extract_seq_dt(cn_obj)

## 4 letters as a group for short, mid, long, long long segments
# A B C D for copy number 0
# E F G H for copy number 1
# I J K L for copy number 2
# M N O P for copy number 3
# Q R S T for copy number 4
# U V W X for copy number 5+

sub_list <- build_sub_matrix()
x$Seqs <- sub_list$map[paste0(x$lenVal, x$segVal)]
x$Seqs

seqs <- extract_sequences(x)

length(seqs$keep)
system.time(
  score_mat <- get_score_matrix(seqs$keep[1:1000], sub_list$mat)
)

system.time(
  score_mat2 <- get_score_matrix(seqs$keep, sub_list$mat, verbose = FALSE)
)

system.time(
  score_mat3 <- get_score_matrix(seqs$keep[1:1000], sub_list$mat, method = "bigmemory")
)

get_score_matrix(seqs$keep[1:10], sub_list$mat, method = "base") %>% typeof()
get_score_matrix(seqs$keep[1:10], sub_list$mat, method = "ff")
get_score_matrix(seqs$keep[1:10], sub_list$mat, method = "bigmemory")[] %>% class() # Return class as base

for (i in LETTERS[1:24]) {
  message("Number of seqs starts with ", i, ": #", sum(startsWith(seqs$keep, i)))
}

set.seed(123)
x <- sapply(1:10000, function(x) {
  paste(sample(LETTERS[1:24], 5, replace = TRUE), collapse = "")
})
y1 <- get_score_matrix(x, sub_list$mat)

y2 <- get_score_matrix2(x, sub_list$mat, verbose = TRUE)

all.equal(y1, y2)

y3 <- get_score_matrix2(x, sub_list$mat, verbose = TRUE, block_size = 3)

z = get_score_matrix2(x, sub_list$mat, cores = 4)
z2 = get_score_matrix2(x, sub_list$mat, cores = 1, verbose = FALSE)
z3 = get_score_matrix(x, sub_list$mat, verbose = FALSE)

all.equal(z2, z3)
all.equal(z, z2)
all.equal(z[, 1:10], z2[, 1:10])
all.equal(z[, 991:1000], z2[, 991:1000])

tr <- sapply(chunk2(1:1000, 100), function(x) {
  all.equal(z[, x], z2[, x])
})

tr

l <- chunk2(1:1000, 100)

q = show_seq_logo(sapply(split(LETTERS[1:24], 1:4), function(x) paste0(x, collapse = "")))
q
show_seq_logo(sapply(split(LETTERS[1:24], 1:4), function(x) paste0(x, collapse = "")))
