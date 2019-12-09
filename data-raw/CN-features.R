## code to prepare `CN.features` dataset goes here
#
# BP10MB: 0、1、2、3、4、5、>5
# BPArm: 0、1、2、3、4、5、6、7、8、9、10、>10<=20、>20<=30、>30
# CN: 0、1、2、3、4、>4<=8、>8
# CNCP: 0、1、2、3、4、>4<=8、>8
# OsCN: 0、1、2、3、4、>4<=10、>10
# (log10 based)SegSize: <=2、>2<=3、>3<=4、>4<=5、>5<=6、>6<=7、>7<=8, >8
# NChrV: 0、1、2、3、4、5、6、7、8、9、10、11、12、13、14、15、16、17、18、19、20、21、22（remove XY）
#
# Length of chains of oscillating copy-number (OsCN)
# Copy number change point (CNCP)
# SegSize (SS)

features = c("BP10MB", "BPArm", "CN", "CNCP", "OsCN", "SS", "NChrV")

# When component is a unique value, min=max
# When component is a data range, make left open and right closed
CN.features = dplyr::tibble(
  feature = features[c(rep(1, 7), rep(2, 14), rep(3, 7), rep(4, 7), rep(5, 7), rep(6, 8), rep(7, 23))],
  min = c(c(0:5, 5),
          c(0:10, 10, 20, 30),
          c(0:4, 4, 8),
          c(0:4, 4, 8),
          c(0:4, 4, 10),
          c(-Inf, 2:8),
          c(0:22)),
  max = c(c(0:5, Inf),
          c(0:10, 20, 30, Inf),
          c(0:4, 8, Inf),
          c(0:4, 8, Inf),
          c(0:4, 10, Inf),
          c(2:8, Inf),
          c(0:22))
)


CN.features = get_feature_components(CN.features)

usethis::use_data(CN.features)