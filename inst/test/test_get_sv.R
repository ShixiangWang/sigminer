load("data-raw/sv.RData")
sv_dt <- read_sv(input = sv)
# usethis::use_package("copynumber", "Suggests")
sv_sig_tally <- sig_tally.SV(sv_dt)
