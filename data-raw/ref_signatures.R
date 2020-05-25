library(tidyverse)
## This has been stored in Maftools
## SBS = vroom::vroom("data-raw/SigProfiler Signatures/SigProfiler reference signatures/SigProfiler reference whole-genome signatures/sigProfiler_SBS_signatures_2018_03_28.csv")
##

## DBS
DBS <- readr::read_csv("data-raw/SigProfiler Signatures/SigProfiler reference signatures/SigProfiler reference whole-genome signatures/sigProfiler_DBS_signatures.csv")
##
DBS_db <- list()
DBS <- DBS %>% column_to_rownames("Mutation Type")
## Ref link: https://cancer.sanger.ac.uk/cosmic/signatures/DBS
DBS_aetiology <- data.frame(
  aetiology = c(
    "UV exposure",
    "exposure to tobacco (smoking) mutagens",
    "Polymerase epsilon exonuclease domain mutations",
    "Unknown",
    "Prior chemotherapy treatment with platinum drugs",
    "Unknown",
    "Defective DNA mismatch repair",
    "Unknown",
    "Unknown",
    "Defective DNA mismatch repair",
    "Possibly related to APOBEC mutagenesis"
  )
)

rownames(DBS_aetiology) <- colnames(DBS)

DBS <- apply(DBS, 2, function(x) x / sum(x)) %>% as.data.frame()

DBS_db <- list(
  db = DBS,
  aetiology = DBS_aetiology,
  date = "2020/05/03"
)

saveRDS(DBS_db, file = "inst/extdata/DBS_signatures.rds")

## INDEL
ID <- readr::read_csv("data-raw/SigProfiler Signatures/SigProfiler reference signatures/SigProfiler reference whole-genome signatures/sigProfiler_ID_signatures.csv")
ID <- ID %>% column_to_rownames("Mutation Type")
rownames(ID)

map_dt <- dplyr::tibble(
  types = rownames(ID)
)

new_rownames <- map_dt %>%
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

rownames(ID) <- new_rownames

## Ref link: https://cancer.sanger.ac.uk/cosmic/signatures/ID
ID_aetiology <- data.frame(
  aetiology = c(
    "Slippage during DNA replication of the replicated DNA strand",
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
    "Unknown"
  )
)

rownames(ID_aetiology) <- colnames(ID)
ID <- apply(ID, 2, function(x) x / sum(x)) %>% as.data.frame()

ID_db <- list(
  db = ID,
  aetiology = ID_aetiology,
  date = "2020/05/03"
)

saveRDS(ID_db, file = "inst/extdata/ID_signatures.rds")

## TSB: Transcriptional Strand Bias Signatures
TSB <- readr::read_csv("data-raw/SigProfiler Signatures/SigProfiler reference signatures/Sigprofiler Transcriptional Strand Bias Signatures/sigProfiler_TSB_signatures.csv")
TSB <- TSB %>%
  mutate(component = paste(Strand,
                           paste0(
                             substr(Subtype, 1, 1),
                             "[",
                             Type,
                             "]",
                             substr(Subtype, 3, 3)
                           ), sep = ":"))

TSB_db <- list()
TSB <- TSB %>%
  dplyr::select(-c("Strand", "Type", "Subtype")) %>%
  column_to_rownames("component")
## Ref link: https://cancer.sanger.ac.uk/cosmic/signatures/SBS
sbs_file = system.file("extdata", "SBS_signatures.RDs",
                       package = "maftools", mustWork = TRUE
)
sbs = readRDS(sbs_file)
rbind(colnames(sbs$db), colnames(TSB))

colnames(TSB) = colnames(TSB) %>% str_remove("-E")

TSB_aetiology <- sbs$aetiology[!rownames(sbs$aetiology) %in% c("SBS84", "SBS85"), , drop = FALSE]

TSB <- apply(TSB, 2, function(x) x / sum(x)) %>% as.data.frame()

TSB_db <- list(
  db = TSB,
  aetiology = TSB_aetiology,
  date = "2020/05/25"
)

saveRDS(TSB_db, file = "inst/extdata/TSB_signatures.rds")
