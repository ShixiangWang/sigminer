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

### Test and Validation

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)
data = data.table::fread("../MatrixGenerator-Test/data/wt/luad/input/test_luad.maf")
colnames(data) = colnames(laml@data)
colnames(data)[14:17] = c("X1", "X2", "Tumor_Sample_Barcode", "X3")

test = read_maf(data)
test_tally = sig_tally(
  test,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  use_syn = TRUE,
  add_trans_bias = TRUE,
  mode = "DBS"
)

test_mat = test_tally$all_matrice$DBS_186

sapply(test_tally$all_matrices, sum)

val_dt = data.table::fread("../MatrixGenerator-Test/data/wt/luad/output/DBS/test.DBS186.all")
val_mat = val_dt %>%
  tibble::column_to_rownames("MutationType") %>%
  as.matrix() %>%
  t

val_mat

dim(test_mat)
dim(val_mat)

val_mat = val_mat[rownames(test_mat), colnames(test_mat)]
all.equal(test_mat, val_mat)

j = 1
for (i in 1:500) {
  dif = sum(abs(test_mat[i, ] -  val_mat[i, ]))
  if (dif > 0) {
    print(dif)
    message(i, " # ", j)
    j = j + 1
  }
}


# Use TCGA-05-4244-01A-01D-1105-08 for deeper debugging
test_mat["TCGA-05-4244-01A-01D-1105-08",]
val_mat["TCGA-05-4244-01A-01D-1105-08",]

rbind(test_mat[60, ],
      val_mat[60, ])

one = read_maf(data[Tumor_Sample_Barcode == "TCGA-05-4244-01A-01D-1105-08"])
one_tally = sig_tally(
  one,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  use_syn = TRUE,
  mode = "ID"
)

sum(one_tally$all_matrice$ID_83)
sum(val_mat["TCGA-05-4244-01A-01D-1105-08",])

rbind(
  one_tally$all_matrice$ID_28,
  val_mat["TCGA-05-4244-01A-01D-1105-08", colnames(one_tally$all_matrice$ID_28), drop = FALSE]
)
