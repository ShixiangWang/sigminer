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
