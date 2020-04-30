library(tidyverse)
## This has been stored in Maftools
## SBS = vroom::vroom("data-raw/SigProfiler Signatures/SigProfiler reference signatures/SigProfiler reference whole-genome signatures/sigProfiler_SBS_signatures_2018_03_28.csv")
##
DBS = vroom::vroom("data-raw/SigProfiler Signatures/SigProfiler reference signatures/SigProfiler reference whole-genome signatures/sigProfiler_DBS_signatures.csv")
##
DBS_db = list()
DBS = DBS %>% column_to_rownames("Mutation Type")
## Ref link: https://cancer.sanger.ac.uk/cosmic/signatures/DBS
DBS_aetiology = data.frame(
  aetiology = c("UV exposure",
                "exposure to tobacco (smoking) mutagens",
                "Polymerase epsilon exonuclease domain mutations",
                "Unknown",
                "Prior chemotherapy treatment with platinum drugs",
                "Unknown",
                "Defective DNA mismatch repair",
                "Unknown",
                "Unknown",
                "Defective DNA mismatch repair",
                "Possibly related to APOBEC mutagenesis")
)

DBS_db = list(
  db = DBS,
  aetiology = DBS_aetiology,
  date = "2020/04/30"
)

saveRDS(DBS_db, file = "inst/extdata/DBS_signatures.rds")

## INDEL
ID = vroom::vroom("data-raw/SigProfiler Signatures/SigProfiler reference signatures/SigProfiler reference whole-genome signatures/sigProfiler_ID_signatures.csv")
ID = ID %>% column_to_rownames("Mutation Type")
rownames(ID)

map_dt = dplyr::tibble(
  types = rownames(ID)
)

new_rownames = map_dt %>%
  mutate(
    types2 = str_remove_all(types, "\\+"),
    is_del = substr(types2, 1, 3) == "DEL",
    is_repeats = str_detect(types2, "_repeats_"),
    is_mh = str_detect(types2, "_MH_"),
    counts = substr(types2, nchar(types2) - 2, nchar(types2))
  ) %>%
  tidyr::separate(counts, into = c("n_ID", "n_size")) %>%
  mutate(
    ID = ifelse(is_del, "Del", "Ins"),
    label = case_when(
      is_mh ~ "M",
      is_repeats ~ "R",
      TRUE ~ substr(types2, 5, 5)
    ),
    new_types = paste(n_ID, ID, label, n_size, sep = ":")
  ) %>%
  pull(new_types)

rownames(ID) = new_rownames

## Ref link: https://cancer.sanger.ac.uk/cosmic/signatures/ID
ID_aetiology = data.frame(
  aetiology = c("Slippage during DNA replication of the replicated DNA strand",
                "Slippage during DNA replication of the template DNA strand",
                "Associated with tobacco smoking",
                "Unknown",
                "Unknown",
                "Defects in DNA-DSB repair by HR",
                "Defective DNA mismatch repair",
                "Unknown. DNA-DSB repair by NHEJ mechanisms",
                "Unknown",
                "Unknown",
                "Unknown",
                "Unknown",
                "UV exposure",
                "Unknown",
                "Unknown",
                "Unknown",
                "Unknown")
)

ID_db = list(
  db = ID,
  aetiology = ID_aetiology,
  date = "2020/04/30"
)

saveRDS(ID_db, file = "inst/extdata/ID_signatures.rds")
