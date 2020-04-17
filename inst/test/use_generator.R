# # install
# # doc: https://github.com/AlexandrovLab/SigProfilerMatrixGeneratorR
# #
# #pip install SigProfilerMatrixGenerator
# #remotes::install_github("AlexandrovLab/SigProfilerMatrixGeneratorR")
#
# library("reticulate")
#
# use_python("/Users/wsx/anaconda3/bin/python", required = TRUE)
# py_config()
#
# #import("SigProfilerMatrixGenerator")
#
# library("SigProfilerMatrixGeneratorR")
# install('GRCh37', rsync=FALSE, bash=TRUE)


# Usage
library("SigProfilerMatrixGeneratorR")
matrices <- SigProfilerMatrixGeneratorR("BRCA", "GRCh37", "/Users/ebergstr/Desktop/BRCA/",
                                        plot=T, exome=F, bed_file=NULL,
                                        chrom_based=F, tsb_stat=False, seqInfo=False, cushion=100)

###

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)
if (require("BSgenome.Hsapiens.UCSC.hg19")) {
  mt_tally <- sig_tally(
    laml,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
    use_syn = TRUE, ignore_chrs = "23"
  )
} else {
  message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
}

val = maftools::trinucleotideMatrix(laml,
                                    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                                    prefix = "chr", add = TRUE,
                                        useSyn = TRUE)


all.equal(mt_tally$nmf_matrix, val$nmf_matrix)
dim(mt_tally$nmf_matrix)
dim(val$nmf_matrix)

all.equal(mt_tally$APOBEC_scores, val$APOBEC_scores)
colnames(mt_tally$APOBEC_scores)
colnames(val$APOBEC_scores)

colnames(mt_tally$APOBEC_scores)[!colnames(mt_tally$APOBEC_scores) == colnames(val$APOBEC_scores)]
colnames(val$APOBEC_scores)[!colnames(mt_tally$APOBEC_scores) == colnames(val$APOBEC_scores)]

purrr::map2(mt_tally$APOBEC_scores, val$APOBEC_scores, all.equal)

all(colnames(mt_tally$nmf_matrix) == colnames(val$nmf_matrix))
z = mt_tally$nmf_matrix[rownames(val$nmf_matrix), colnames(val$nmf_matrix)]

all.equal(z, val$nmf_matrix)
setdiff(rownames(mt_tally$nmf_matrix), rownames(val$nmf_matrix))

"TCGA-AB-2823" %in% laml@data$Tumor_Sample_Barcode
"TCGA-AB-2823" %in% rownames(val$nmf_matrix)

## Test 1536
data = data.table::fread("../MatrixGenerator-Test/data/wt/luad/input/test_luad.maf")
colnames(data) = colnames(laml@data)
colnames(data)[14:17] = c("X1", "X2", "Tumor_Sample_Barcode", "X3")

test = read_maf(data)
test_tally = sig_tally(
  test,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  use_syn = TRUE
)

test_mat = test_tally$all_matrice$SBS_1536

table(data$Chromosome)

val_dt = data.table::fread("../MatrixGenerator-Test/data/wt/luad/output/SBS/test.SBS1536.all")
val_mat = val_dt %>%
  tibble::column_to_rownames("MutationType") %>%
  as.matrix() %>%
  t

val_mat

dim(test_mat)
dim(val_mat)

all(colnames(test_mat) == colnames(val_mat))

val_mat = val_mat[rownames(test_mat), colnames(test_mat)]

all.equal(test_mat, val_mat)

for (i in 1:nrow(test_mat)) {
  if (sum(test_mat[i, ] -  val_mat[i, ]) > 0) {
    print(i)
  }
}

sum(test_mat)
sum(val_mat)
data[Variant_Type == "SNP"] %>% nrow()

which(test_mat[4, ] - val_mat[4, ] > 0)
test_mat[4, 183]
val_mat[4, 183]

sum(test_tally$all_matrice$SBS_6)
sum(test_tally$all_matrice$SBS_96)
