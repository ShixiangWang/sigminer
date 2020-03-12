# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))

# Prepare copy number signature analysis
options(sigminer.sex = "male", sigminer.copynumber.max = 10)

cn_prepare <- sig_tally(cn, method = "W")
cn_prepare <- sig_tally(cn, method = "W", cores = 2)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89),
  segVal = c(2, 2, 1, 2, 3, 3, 2, 2, 2)
)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89),
  segVal = c(2, 2, 2, 2, 3, 3, 2, 2, 2)
)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89),
  segVal = c(1, 2, 1, 2, 3, 3, 2, 2, 2)
)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89, 99),
  segVal = c(2, 2, 1, 2, 3, 3, 2, 2, 2, 1)
)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89, 99),
  segVal = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89, 99),
  segVal = c(1, 2, 2, 1, 2, 1, 2, 3, 2, 2)
)


df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60),
  end = c(9, 19, 29, 39, 49, 59, 69),
  segVal = c(2, 1, 3, 3, 3, 2, 2)
)

df

join_segments(df)
join_segments(df[c(1:2), ])
join_segments(df[c(1), ])
join_segments(df[c(2:4), ])
join_segments(df[c(5:8), ])


library(maftools)
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml = read.maf(maf = laml.maf)

library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
laml.tnm = trinucleotideMatrix(maf = laml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
dim(laml.tnm$nmf_matrix)
library('NMF')
laml.sig = extractSignatures(mat = laml.tnm, n = 6, pConstant = 1e-9)
laml.se = signatureEnrichment(maf = laml, sig_res = laml.sig)
sum(laml.se$cf_sizes$N)

load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
# Extract copy number signatures
library(NMF)
sig <- sig_extract(cn_prepare$nmf_matrix, 5, nrun = 10)

# Methods 'consensus' and 'samples' are from NMF::predict()
get_groups(sig)
get_groups(sig, method = "exposure")
get_groups(sig, method = "consensus")
get_groups(sig, method = "consensus", match_consensus = TRUE)
get_groups(sig, method = "samples")

debug(get_groups)
undebug(get_groups)

## Test New feature

# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))
# Use method designed by Wang, Shixiang et al.
cn_prepare <- sig_tally(cn, method = "W", feature_setting = CN.features[1:50])
cn_prepare <- sig_tally(cn, method = "W", feature_setting = CN.features)

cn_prepare$features$BoChr
cn_prepare$features$NChrV
cn_prepare$features$NC50
cn_prepare$components
cn_prepare$parameters
debug(sig_tally)

library(NMF)
sigs <- sig_extract(cn_prepare$nmf_matrix, n_sig = 3, pConstant = 1e-9)
show_sig_profile(sigs, method = "W", normalize = "feature", style = "cosmic")
show_sig_profile(sigs, method = "W", normalize = "feature", style = "cosmic", filters = c("NChrV", "BoChr"))
show_sig_profile(sigs, method = "W", normalize = "feature", style = "cosmic", filters = c("BoChr"))

# Load copy number prepare object
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
# Extract copy number signatures
library(NMF)
res <- sig_extract(cn_prepare$nmf_matrix, 5, nrun = 2)
show_sig_profile(res, normalize = "feature", style = "cosmic")
show_sig_profile(res, normalize = "feature", style = "cosmic", filters = c("bp10MB", "bpchrarm"))
show_sig_profile(laml.sig$signatures, mode = "SBS", style = "default",rm_panel_border = F)
show_sig_profile(laml.sig$signatures, mode = "SBS", style = "cosmic")

p = ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
add_h_arrow(p, x = 3, y = 30)



# Change name of function
#sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl R/*`
#sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl docs/*`
#sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl man/*`
#sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl tests/*`
#sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl _pkgdown.yml`


## Test copy number exposure
# Load toy dataset of absolute copynumber profile
load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))
cn <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE, verbose = TRUE
)
cn.tally <- sig_tally(cn, method = "W", feature_setting = CN.features)
library(NMF)
cn.sig <- sig_extract(cn.tally$nmf_matrix, n_sig = 2, pConstant = 1e-9)
sig.expo <- get_sig_exposure(cn.sig)

sig.expo[, Total := Sig1 + Sig2]

sig.expo
cn_dt = cn@data[, .(total = .N), by = sample]
dt = merge(cn_dt, sig.expo, by = "sample")

cor(dt$total, dt$Total, method = "spearman")
plot(dt$total, dt$Total, xlim = c(25, 55), ylim = c(25, 55))


## An sample for mutational signature
library(maftools)
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml = read.maf(maf = laml.maf)

library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
laml.tnm = trinucleotideMatrix(maf = laml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
dim(laml.tnm$nmf_matrix)

library(NMF)
library(sigminer)

laml.sig = sig_extract(laml.tnm$nmf_matrix, n_sig = 3, pConstant = 1e-9)
get_sig_exposure(laml.sig)
get_sig_exposure(laml.sig, type = "relative")


## Group exposure
# Load copy number signature
load(system.file("extdata", "toy_copynumber_signature.RData",
                 package = "sigminer", mustWork = TRUE
))
# Show signature exposure
show_sig_exposure(sig, groups = c(rep("grp1", 5), rep("grp2", 5)))
show_sig_exposure(sig, groups = c(rep("grp1", 5), rep("grp2", 5)), style = "cosmic")
show_sig_exposure(sig, groups = c(rep("grp1", 5), rep("grp2", 5)), grp_order = c("grp1", "grp2"))
show_sig_exposure(sig, groups = c(rep("grp1", 5), rep("grp2", 5)), grp_order = c("grp2", "grp1"), style = "cosmic")
show_sig_exposure(sig, groups = c(rep("grp1", 5), rep("grp2", 5)), grp_order = c("grp2", "grp1"), style = "cosmic",
                  legend_position = "none")
show_sig_exposure(sig, groups = c(rep("grp1", 5), rep("grp2", 5)), grp_order = c("grp2", "grp1"), style = "cosmic",
                  legend_position = "none", rm_space = TRUE, rm_grid_line = FALSE, rm_panel_border = T)


## Circos for copy number profile
library(circlize)

load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))
cn <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE, verbose = TRUE
)


# set.seed(999)
# bed = generateRandomBed(nr = 100)
#
# circos.initializeWithIdeogram()
# circos.genomicTrack(bed, panel.fun = function(region, value, ...) {
#   print(region)
#   print(value)
#   circos.genomicPoints(value, col = ifelse(value[[1]] > 0, "red", "green"))
# })
#
#
# set.seed(123)
# bed1 = generateRandomBed(nr = 100)
# bed1 = bed1[sample(nrow(bed1), 20), ]
# bed2 = generateRandomBed(nr = 100)
# bed2 = bed2[sample(nrow(bed2), 20), ]
#
# circos.initializeWithIdeogram()
# circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5),
#                    border = NA)

segTabs$sample = NULL
segTabs$chr = paste0("chr", segTabs$chr)
colnames(segTabs)[1] = "chr"
circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", c(3,5,2,8)))
circos.initializeWithIdeogram(species = "hg19")
# bed = generateRandomBed(nr = 100, nc = 1)
col_fun = colorRamp2(c(1, 2, 4), c("blue", "black", "red"))
circos.genomicHeatmap(segTabs, col = col_fun, side = "inside", border = "white", border_lwd = 0)
circos.clear()


## Scoring
# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))



devtools::install_github("ShixiangWang/sigminer@V3")

library(sigminer)
library(NMF)

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)
if (require("BSgenome.Hsapiens.UCSC.hg19")) {
  mt_prepare <- sig_tally(
    laml,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
    prefix = "chr", add = TRUE, useSyn = TRUE
  )
} else {
  message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
}

sig = sig_extract(mt_prepare$nmf_matrix, n_sig = 3, pConstant = 1e-13)
get_sig_exposure(sig)
