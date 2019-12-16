library(maftools)
?read.maf
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read.maf(maf = laml.maf)
laml.tnm <- trinucleotideMatrix(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
                                prefix = 'chr', add = TRUE, useSyn = TRUE)
laml.sign <- extractSignatures(mat = laml.tnm, n = 3, pConstant = 0.001)
laml.sign$signatures

compareSignatures(laml.sign)


load(system.file("extdata", "toy_copynumber_signature.RData",
                 package = "sigminer", mustWork = TRUE
))
# Load copy number prepare object
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
params = get_tidy_parameter(cn_prepare$components)
show_sig_profile(sig, params = params, y_expand = 2)
show_sig_profile(sig, params = params, y_expand = 3, show_cv = T, params_label_angle = 70)



# Optimize signature exposure
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
res <- sig_auto_extract(cn_prepare2$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1)

str(cn_prepare$features$segsize)


load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))

segTabs$start = as.character(segTabs$start)
segTabs$end = as.numeric(segTabs$end)

system.time(
  cn <- read_copynumber(segTabs,
                        seg_cols = c("chromosome", "start", "end", "segVal"),
                        genome_build = "hg19", complement = FALSE, verbose = TRUE
  )
)


# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))
# Prepare copy number signature analysis
system.time(
  cn_prepare3 <- sig_derive(cn, method = "M")
)
system.time(
  cn_prepare2 <- sig_derive(cn, method = "W")
)
system.time(
  cn_prepare3 <- sig_derive(cn, method = "W")
)

str(cn_prepare$features$copynumber)

cc = getBPnum(cn_list, chrlen = sigminer::chromsize.hg19)
cc = getOscilation(cn_list)
cc2 = getOscilation(cn_list)

cc = getCentromereDistCounts(cn_list, centromeres.hg19, chromsize.hg19)
cc2 = getCentromereDistCounts2(cn_list, centromeres.hg19, chromsize.hg19)
debug(getCentromereDistCounts2)


show_cn_features(cn_prepare$features)
show_cn_features(cn_prepare2$features[1:6], method = "Wang")

show_cn_components(cn_prepare$parameters)
show_cn_components(cn_prepare2$parameters, method = "W")

# Test show signature profile
sigs = sig_auto_extract(cn_prepare2$nmf_matrix)

library(NMF)
sigs_est = sig_estimate(cn_prepare2$nmf_matrix, pConstant = 0.001, verbose = TRUE)
show_rank_survey(sigs_est)
sigs = sig_extract(cn_prepare2$nmf_matrix, n_sig = 3, pConstant = 0.001)

sigs2 = sig_extract(cn_prepare$nmf_matrix, n_sig = 3)

show_sig_profile(sigs, method = "W")
show_sig_profile(sigs2)

sig = sig_auto_extract(cn_prepare2$nmf_matrix[, 1:50])
show_sig_profile(sig, method = "W")

# Load toy dataset of absolute copynumber profile
load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))
cn <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE, verbose = TRUE
)

cc = sig_derive(cn, method = "W")
ss = sig_auto_extract(cc$nmf_matrix)
show_sig_profile(ss, method = "W")
show_sig_exposure(ss)
get_sig_exposure(ss, type = "relative")

show_cn_features(cc$features, method = "W")
show_cn_components(cc$parameters, method = "W")


data <- dplyr::tibble(
  Group1 = rep(LETTERS[1:5], each = 10),
  Group2 = rep(LETTERS[6:15], each = 5),
  zzzz = c(rep("xx", 20), rep("yy", 20), rep(NA, 10))
)

show_group_mapping(data, sig_col = "Group1", map_cols = colnames(data)[-1], fill_na = "NA")

show_group_mapping(data, sig_col = "Group1", map_cols = colnames(data)[-1], include_sig = TRUE, fill_na = "NA")

df %>% mutate_all(funs(ifelse(is.na(.), replace_na(., 0), .)))


# Load copy number signature
load(system.file("extdata", "toy_copynumber_signature.RData",
                 package = "sigminer", mustWork = TRUE
))
# Show signature profile
show_sig_profile(sig, normalize = "column")


c(
  "cyan", "red", "yellow", "purple",
  "green", "blue", "black", "gray"
)
#1 GnBu
#2 Reds
#3 Oranges
#4 Purples
#5 Greens
#6 Blues
#7 Greys
colorRampPalette(RColorBrewer::brewer.pal(3, name = "Blues"))(50)

# Association analysis
library(ggcor)
library(vegan) # 使用vegan包所带的数据集
data(varechem)
data(varespec)
corr <- fortify_cor(varechem, type = "upper", show.diag = TRUE,
                    cor.test = TRUE, cluster.type = "all")

mantel <- fortify_mantel(varespec[, 1:4], varechem,
                         spec.select = list(spec01 = 1,
                                            spec02 = 2,
                                            spec03 = 3,
                                            spec04 = 4),
                         mantel.fun = "mantel.randtest")
group <- rep(LETTERS[1:3], 8)
mantel = fortify_mantel(varespec[ , 38:43], varechem,
                        spec.group = group,
                        env.group = group,
                        is.pair = TRUE,
                        mantel.fun = "mantel.randtest")

ggcor(corr, xlim = c(-5, 14.5)) +
  add_link(mantel, diag.label = TRUE) +
  add_diaglab(angle = 45) +
  geom_square() + remove_axis("y")


corr <- cor(mtcars)
df <- as_cor_tbl(corr)
?as_cor_tbl_fct()
