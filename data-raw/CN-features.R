## code to prepare `CN.features` dataset goes here
# feature     components   dist      mean      sd n_obs
# <chr>       <chr>        <chr>    <dbl>   <dbl> <dbl>
#   1 bp10MB      bp10MB1      pois  3.66e-12 1.91e-6  2839
# 2 bp10MB      bp10MB2      pois  9.24e- 1 9.61e-1   161
# 3 bpchrarm    bpchrarm1    pois  1.22e- 6 1.11e-3   328
# 4 bpchrarm    bpchrarm2    pois  1.21e+ 0 1.10e+0   102
# 5 bpchrarm    bpchrarm3    pois  5.54e+ 0 2.35e+0    10
# 6 changepoint changepoint1 norm  9.01e- 1 7.10e-1   236
#
# BP10MB：0、1、2、3、4、5、>5
# CopyNumber：0、1、2、3、4、>4<=8、>8
# ChangePoint： 0、1、2、3、4、>8
# BPChrArm：0、1、2、3、4、5、6、7、8、9、10、>10<=20、>20<=30、>30
# OSCN：0、1、2、3、4、>4<=10、>10
# Log10SegSize（log10）：<=2、>2<=3、>3<=4、>4<=5、>5<=6、>6<=7、>7<=8
# NChrV：0、1、2、3、4、5、6、7、8、9、10、11、12、13、14、15、16、17、18、19、20、21、22（0 to maximun chr number）
# Length of segments with oscillating copy-number
# Change in copy number

features = c("BP10MB", "BPChrArm", "CopyNumber", "ChangePoint", "LenOsCN", "Log10SegSize", "NChrV")

usethis::use_data("CN.features")
