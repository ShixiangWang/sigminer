## Signature identification for absolute copy number profile

load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))
cn <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE, verbose = TRUE
)

save(cn, file = "inst/extdata/toy_copynumber.RData")

# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))

# Prepare copy number signature analysis
options(sigminer.sex = "male", sigminer.copynumber.max = 10)

cn_tally_M = sig_tally(cn, method = "M")
cn_tally_W = sig_tally(cn, method = "W")

save(cn_tally_M, file = "inst/extdata/toy_copynumber_tally_M.RData")
save(cn_tally_W, file = "inst/extdata/toy_copynumber_tally_W.RData")


load(file = "inst/extdata/toy_copynumber_tally_W.RData")
library(NMF)
sig = sig_extract(cn_tally_W$nmf_matrix, n_sig = 2, pConstant = 1e-13)
save(sig, file = "inst/extdata/toy_copynumber_signature_by_W.RData")


sig2 = sig_extract(cn_tally_W$nmf_matrix, n_sig = 2, pConstant = 1e-13, optimize = TRUE)

identical(sig$Signature, sig2$Signature)
identical(sig$Signature.norm, sig2$Signature.norm)

sum((apply(sig$Signature[22:28,], 2, function(x) x/sum(x)) %*% sig$Exposure - t(cn_tally_W$nmf_matrix)[22:28, ])^2)
sum((apply(sig2$Signature[22:28,], 2, function(x) x/sum(x)) %*% sig2$Exposure - t(cn_tally_W$nmf_matrix)[22:28, ])^2)

expo1 = get_sig_exposure(sig)
expo2 = get_sig_exposure(sig2)

plot(expo1$Sig1, expo2$Sig1)
plot(expo1$Sig2, expo2$Sig2)

cor(expo1$Sig1, expo2$Sig1)
cor(expo1$Sig2, expo2$Sig2)

## Signature identification for SBS
library(sigminer)
library(NMF)

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)
laml <- read_maf(maf = laml.maf)
if (require("BSgenome.Hsapiens.UCSC.hg19")) {
  mt_tally <- sig_tally(
    laml,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
    prefix = "chr", add = TRUE, useSyn = TRUE
  )
} else {
  message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
}

sig = sig_extract(mt_tally$nmf_matrix, n_sig = 3, pConstant = 1e-13)
get_sig_exposure(sig)


sig2 = sig_extract(mt_tally$nmf_matrix, n_sig = 3, pConstant = 1e-13, optimize = TRUE)

identical(sig$Signature, sig2$Signature)
identical(sig$Signature.norm, sig2$Signature.norm)

sum((apply(sig$Signature, 2, function(x) x/sum(x)) %*% sig$Exposure - t(mt_tally$nmf_matrix))^2)
sum((apply(sig2$Signature, 2, function(x) x/sum(x)) %*% sig2$Exposure - t(mt_tally$nmf_matrix))^2)

expo1 = get_sig_exposure(sig)
expo2 = get_sig_exposure(sig2)

plot(expo1$Sig1, expo2$Sig1)
plot(expo1$Sig2, expo2$Sig2)

cor(expo1$Sig1, expo2$Sig1)
cor(expo1$Sig2, expo2$Sig2)

show_sig_profile(mt_tally$nmf_matrix[1, , drop = FALSE] %>% t(), mode = "SBS",
                 style = "cosmic", normalize = "raw",
                 check_sig_names = FALSE)

show_sig_profile(mt_tally$nmf_matrix[1:2, ] %>% t(), mode = "SBS",
                 style = "cosmic", normalize = "raw",
                 check_sig_names = FALSE)

show_cosmic_sig_profile(sig_index = 1, style = "cosmic")


# Change name of function
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl R/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl docs/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl man/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl tests/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl _pkgdown.yml`


load(system.file("extdata", "asso_data.RData",
                 package = "sigminer", mustWork = TRUE
))

#?show_sig_feature_corrplot
show_sig_feature_corrplot(tidy_data.seqz.feature, p_val = 1)
show_sig_feature_corrplot(tidy_data.seqz.feature, breaks_count = NULL)
show_sig_feature_corrplot(tidy_data.seqz.feature, breaks_count = NA)


# Test DBS ----------------------------------------------------------------

maf = read_maf("../MatrixGenerator-Test/data/LUSC/tcga_lusc_from_dcc.maf")
mt_tally <- sig_tally(
  maf,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  use_syn = TRUE, add_trans_bias = TRUE, mode = "DBS"
)



z2 = query[, search_DBS(.SD),
      by = Tumor_Sample_Barcode]

# query %>%
#   dplyr::as_tibble() %>%
#   dplyr::group_by(Tumor_Sample_Barcode) %>%
#   tidyr::nest() %>%
#   dplyr::mutate(data = purrr::map(.data$data, search_DBS)) %>%
#   tidyr::unnest("data") %>%
#   data.table::as.data.table()


# Test INDEL --------------------------------------------------------------
debug(generate_matrix_INDEL)
mt_tally <- sig_tally(
  maf,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  use_syn = TRUE, mode = "ID", add_trans_bias = TRUE
)


# Test bootstrap batch in practice ----------------------------------------

library(sigminer)
library(NMF)

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)
laml <- read_maf(maf = laml.maf)
if (require("BSgenome.Hsapiens.UCSC.hg19")) {
  mt_tally <- sig_tally(
    laml,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
    prefix = "chr", useSyn = TRUE
  )
} else {
  message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
}

sig_fit(mt_tally$nmf_matrix %>% t(),
                        sig_index = 1:30,
                        sig_db = "legacy",
                        methods = c("LS", "QP"))

for (i in 1:nrow(mt_tally$nmf_matrix)) {
  print(i)
  bt_res <- sig_fit_bootstrap(mt_tally$nmf_matrix[i, ],
                              sig_index = 1:30,
                              sig_db = "legacy",
                              methods = c("QP"), n = 1)
}

bt_res <- sig_fit_bootstrap(mt_tally$nmf_matrix[1, ],
                            sig_index = 1:30,
                            sig_db = "legacy",
                            methods = c("QP"), n = 2)

bt_res <- sig_fit_bootstrap_batch(mt_tally$nmf_matrix %>% t(),
                                  sig_index = 1:30,
                                  sig_db = "legacy",
                                  methods = c("QP"), n = 2)

load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))
debug(get_LengthFraction)
cn <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE
)


cn2 <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE
)


identical(cn@annotation, cn2@annotation)


load("~/biodata/DoAbsolute/CN_list.RData")
cn <- read_copynumber(input = CN_list$ACC@data, genome_build = "hg38",
                      seg_cols = c("chromosome", "start", "end", "segVal"))
cn2 <- read_copynumber(input = CN_list$ACC@data, genome_build = "hg38",
                      seg_cols = c("chromosome", "start", "end", "segVal"))

identical(cn@annotation, cn2@annotation)


# Signature heatmap -------------------------------------------------------

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)

mt_tally <- sig_tally(
  laml,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  use_syn = TRUE, mode = "ALL", add_trans_bias = TRUE
)

str(mt_tally, max.level = 1)

show_sig_profile_heatmap(mt_tally$SBS_96[1:2,] %>% t(), check_sig_names = F)

show_sig_profile_heatmap(mt_tally$SBS_24[1:2,] %>% t(), check_sig_names = F, flip_xy = TRUE)
show_sig_profile_heatmap(mt_tally$SBS_1536[1:2,] %>% t(), check_sig_names = F)
show_sig_profile_heatmap(mt_tally$SBS_384[1:2,] %>% t(), check_sig_names = F, flip_xy = TRUE)
show_sig_profile_heatmap(mt_tally$SBS_6144[1:2,] %>% t(), check_sig_names = F)

debug(show_sig_profile_heatmap)

maf = read_maf("../MatrixGenerator-Test/data/LUSC/tcga_lusc_from_dcc.maf")
mt_tally <- sig_tally(
  maf,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  use_syn = TRUE, add_trans_bias = TRUE, mode = "DBS"
)

show_sig_profile_heatmap(mt_tally$all_matrices$DBS_78[3:4,] %>% t(), check_sig_names = F, mode = "DBS")
show_sig_profile_heatmap(mt_tally$all_matrices$DBS_186[3:4,] %>% t(), check_sig_names = F, mode = "DBS")
show_sig_profile_heatmap(mt_tally$all_matrices$DBS_1248[3:4,] %>% t(), check_sig_names = F, mode = "DBS")


# Test sigprofiler --------------------------------------------------------

load(system.file("extdata", "toy_copynumber_tally_M.RData",
                 package = "sigminer", mustWork = TRUE
))

sigprofiler_extract(cn_tally_M$nmf_matrix, "~/test/test_sigminer", use_conda = TRUE)
sigprofiler_extract(cn_tally_M$nmf_matrix, "~/test/test_sigminer", use_conda = FALSE, py_path = "/Users/wsx/anaconda3/bin/python",
            is_exome = TRUE)
sigprofiler_extract(cn_tally_M$nmf_matrix, "~/test/test_sigminer", use_conda = FALSE, py_path = "/Users/wsx/anaconda3/bin/python",
            is_exome = TRUE, genome_build = "hg38")

sigprofiler_extract(cn_tally_M$nmf_matrix, "~/test/test_sigminer", use_conda = FALSE,
                    py_path = "C:\\Tools\\miniconda3\\envs\\learn\\python.exe",
                    is_exome = TRUE)


# Test enhancement of show_sig_profile ------------------------------------

load(system.file("extdata", "toy_mutational_signature.RData",
                 package = "sigminer", mustWork = TRUE
))
# Show signature profile
p1 <- show_sig_profile(sig2, mode = "SBS", sig_orders = 1)
p2 <- show_sig_profile(sig2, mode = "SBS", sig_orders = 2)

show_sig_profile_loop(sig2, mode = "SBS")
show_sig_profile_loop(sig2, mode = "SBS", style = "cosmic", sig_names = c("A", "B", "C"))


# Test New CN class -------------------------------------------------------

copyNumberInfo = function(chrom,start,end,CN,CNb=NULL,doVariables=TRUE)
{
  # segment lengths
  print("lengths")
  lengths = (end-start)/1000000
  if(!is.null(CNb))
  {
    # allelic imbalance
    print("imbalance")
    CNa = CN-CNb
    imba = abs(CNa-CNb)
    # LOH
    LOH=pmin(CNa,CNb)
  } else {
    imba=NULL
    LOH=NULL
  }
  # combine
  print("combine")
  if(!doVariables)
  {
    combined = list(CN=CN,
                    lengths=lengths,
                    LOH=LOH)
    return(combined)
  }
  LOHstatus = ifelse(LOH==0,"LOH","het")
  LOHstatus[which(CN==0)] = "homdel"
  variables = table(paste(cut(CN,
                              breaks=c(-0.5,1.5,2.5,4.5,Inf),
                              labels=c("del","neut","dup","amp")),
                          LOHstatus,
                          cut(lengths,breaks=c(-0.01,0.01,0.1,1,10,Inf)),
                          sep=":"))
  return(variables)
}

load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))
set.seed(1234)
segTabs$minor_cn <- sample(c(0, 1), size = nrow(segTabs), replace = TRUE)
copyNumberInfo(
  segTabs$chromosome,
  segTabs$start,
  segTabs$end,
  segTabs$segVal,
  CNb=segTabs$minor_cn,
  doVariables=TRUE) -> x

names(x)

x2 <- data.table::fread("inst/extdata/metastatic_tumor.segtab.txt")
copyNumberInfo(
  x2$Chromosome,
  x2$Start.bp,
  x2$End.bp,
  x2$modal_cn,
  x2$subclonal,
  doVariables=TRUE) -> x2

cn <- read_copynumber(
  segTabs,
  seg_cols = c("chromosome", "start", "end", "segVal"),
  genome_build = "hg19", complement = FALSE, add_loh = TRUE,
  join_adj_seg = FALSE)
dt <- cn@data
dt
