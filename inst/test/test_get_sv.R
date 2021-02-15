# sv <- readxl::read_xlsx("inst/extdata/sv.xlsx") %>%
#   as.data.frame()
# save(sv, file = "./inst/extdata/sv.RData")

load(system.file("extdata","sv.RData",package = "sigminer",mustWork = TRUE))
sv_dt <- read_sv(input = sv)
# usethis::use_package("copynumber", "Suggests")
sv_sig_tally <- sig_tally.SV(sv_dt)
