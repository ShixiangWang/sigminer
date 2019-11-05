load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
library(NMF)
cn_estimate <- sig_estimate(cn_prepare$nmf_matrix,
                            cores = 1, nrun = 5,
                            verbose = TRUE
)



## It seems parallel computation makes error on Windows

# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))
# Prepare copy number signature analysis
cn_prepare <- derive(cn)


# Prepare mutational signature analysis
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)
library(BSgenome.Hsapiens.UCSC.hg19)
mt_prepare <- derive(
  laml,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  prefix = "chr", add = TRUE, useSyn = TRUE
)


params = structure(list(feature = c("bp10MB", "bp10MB", "bp10MB", "bpchrarm",
                                    "bpchrarm", "bpchrarm", "bpchrarm", "changepoint", "changepoint",
                                    "changepoint", "changepoint", "copynumber", "copynumber", "copynumber",
                                    "copynumber", "copynumber", "osCN", "osCN", "osCN", "segsize",
                                    "segsize", "segsize", "segsize", "segsize", "segsize", "segsize",
                                    "segsize", "segsize", "segsize", "segsize", "segsize", "segsize",
                                    "segsize", "segsize"), components = c("bp10MB1", "bp10MB2", "bp10MB3",
                                                                          "bpchrarm1", "bpchrarm2", "bpchrarm3", "bpchrarm4", "changepoint1",
                                                                          "changepoint2", "changepoint3", "changepoint4", "copynumber1",
                                                                          "copynumber2", "copynumber3", "copynumber4", "copynumber5", "osCN1",
                                                                          "osCN2", "osCN3", "segsize1", "segsize2", "segsize3", "segsize4",
                                                                          "segsize5", "segsize6", "segsize7", "segsize8", "segsize9", "segsize10",
                                                                          "segsize11", "segsize12", "segsize13", "segsize14", "segsize15"
                                    ), dist = c("pois", "pois", "pois", "pois", "pois", "pois", "pois",
                                                "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm",
                                                "norm", "pois", "pois", "pois", "norm", "norm", "norm", "norm",
                                                "norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm",
                                                "norm", "norm", "norm"), mean = c(5.17530078218012e-09, 0.784423632270936,
                                                                                  3.61203200085818, 5.3948296166391e-12, 1.97434778799569, 7.48484436651666,
                                                                                  20.1864033272021, 1, 1.75343523019778, 2, 8.6779532929398, 0.769290721137994,
                                                                                  2, 3.37814380921201, 5.11770846782, 9.60065618997599, 0.339655037283972,
                                                                                  2.80587520441537, 9.78101019515257, 90852.8920645451, 556896.222888657,
                                                                                  1711028.95829317, 3853161.71330275, 8086167.83548332, 16321712.6291783,
                                                                                  28823161.0035946, 37556698.481876, 40554858.7576789, 45190141.7850124,
                                                                                  62777804.7252656, 81478322.1324386, 89266197.7641066, 133699713.396281,
                                                                                  145648106.44812), sd = c(7.19395633999826e-05, 0.885676934480591,
                                                                                                           1.90053466184076, 2.32267725193129e-06, 1.40511486647736, 2.73584436079918,
                                                                                                           4.49292814623182, 2.22047852496938e-16, 1.92668219056823, 4.44104074014483e-16,
                                                                                                           4.83298714197501, 0.421303619660361, 6.66142673349964e-16, 0.493078206726886,
                                                                                                           1.02326106741676, 4.44005990690756, 0.582799311327641, 1.6750746862201,
                                                                                                           3.12746066244686, 69484.4184390171, 287498.595420422, 718301.143933478,
                                                                                                           1388049.63284035, 2769474.91908646, 4908143.50178269, 5158263.70292865,
                                                                                                           9020577.03474325, 989189.975178436, 10244849.4745123, 7055227.31856478,
                                                                                                           1273751.27630876, 8393648.48018664, 495275.92136265, 47616885.068461
                                                                                  ), n_obs = c(254716, 37886, 2251, 25755, 13454, 2608, 469, 34187,
                                                                                               14700, 14939, 1036, 18847, 37599, 20534, 7744, 1572, 23376, 1293,
                                                                                               124, 4034, 8450, 11643, 10468, 10214, 10870, 4958, 2550, 2259,
                                                                                               3257, 5068, 1907, 3332, 1190, 6096)), class = c("tbl_df", "tbl",
                                                                                                                                               "data.frame"), row.names = c(NA, -34L))
show_cn_components(params, show_weights = F)

params2 = params
params2$sd[params$components == "copynumber2"] = 6.66e-2
show_cn_components(params2, show_weights = T)


#=====================
# Load toy dataset of absolute copynumber profile
load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))
cn <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE, verbose = TRUE
)
