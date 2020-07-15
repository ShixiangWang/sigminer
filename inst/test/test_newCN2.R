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

z[1:10, l[[11]]]
z2[1:10, l[[11]]]

sum(z) - sum(z2)

require(ggplot2)
require(ggseqlogo)

# Some sample data
data(ggseqlogo_sample)

require(ggplot2)
require(ggseqlogo)

# Some sample data
data(ggseqlogo_sample)

p1 = ggseqlogo( seqs_dna$MA0001.1, method = 'bits' )
p1

seqs_numeric = chartr('ATGC','1234', seqs_dna$MA0001.1)
ggseqlogo(seqs_numeric, method='p', namespace=1:4)

ggseqlogo(x[1:10])


cs1 = make_col_scheme(chars=c('A', 'T', 'C', 'G'), groups=c('gr1', 'gr1', 'gr2', 'gr2'),
                      cols=c('purple', 'purple', 'blue', 'blue'))

# Generate sequence logo
ggseqlogo(seqs_dna$MA0001.1, col_scheme=cs1)
