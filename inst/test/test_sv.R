inputdata <- readxl::read_xlsx(path = "F:/sv.xlsx",col_names = T)

sv = read_sv(inputdata)
sv_object <- sig_tally.SV(sv)

svlist <- get_svlist(sv)
getRearrSize_v1(svlist)
test_clustered <- getClustered_v1(svlist)
getType_v1(svlist)
test_feature <- get_features_sv(svlist)
test_component <- get_components_sv(test_feature)
test_matrix <- get_matrix_sv(CN_components = test_component)
