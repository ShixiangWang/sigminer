library(tidyverse)
devtools::load_all()

## SBS
download.file("https://raw.githubusercontent.com/AlexandrovLab/SigProfilerExtractor/master/SigProfilerExtractor/data/sigProfiler_SBS_signatures_2018_03_28.csv",
              "data-raw/sigProfiler_SBS_signatures_2018_03_28.csv")
SBS = readr::read_csv("data-raw/sigProfiler_SBS_signatures_2018_03_28.csv")
SBS = SBS %>%
  dplyr::mutate(
    Type = paste0(
      substr(SubType, 1, 1),
      "[",
      Type,
      "]",
      substr(SubType, 3, 3)
    )
  ) %>%
  dplyr::select(-SubType) %>%
  tibble::column_to_rownames("Type")

cosmic_v3.1 <- readxl::read_excel("data-raw/COSMIC_Mutational_Signatures_v3.1.xlsx", sheet = 1)
cosmic_v3.1 = cosmic_v3.1 %>%
  dplyr::mutate(
    Type = paste0(
      substr(Subtype, 1, 1),
      "[",
      Type,
      "]",
      substr(Subtype, 3, 3)
    )
  ) %>%
  dplyr::select(-Subtype) %>%
  tibble::column_to_rownames("Type")

get_sig_similarity(cosmic_v3.1, SBS)
colSums(cosmic_v3.1)

# sbs = get_sig_db("SBS")
# get_sig_similarity(as.matrix(SBS), sbs$db)
# get_sig_similarity(SBS, sbs$db)
# get_sig_similarity(SBS, as.data.frame(sbs$db))

SBS = cosmic_v3.1

SBS_aetiology <- data.frame(
  aetiology = c(
    "Spontaneous deamination of 5-methylcytosine (clock-like signature)",
    "Activity of APOBEC family of cytidine deaminases",
    "Defective homologous recombination DNA damage repair",
    "Tobacco smoking",
    "Unknown (clock-like signature)",
    "Defective DNA mismatch repair",
    "Ultraviolet light exposure",
    "Ultraviolet light exposure",
    "Ultraviolet light exposure",
    "Ultraviolet light exposure",
    "Unknown",
    "Polimerase eta somatic hypermutation activity",
    "Polymerase epsilon exonuclease domain mutations",
    "Polymerase epsilon exonuclease domain mutations",
    "Temozolomide treatment",
    "Unknown",
    "Activity of APOBEC family of cytidine deaminases",
    "Concurrent polymerase epsilon mutation and defective DNA mismatch repair",
    "Defective DNA mismatch repair",
    "Unknown",
    "Unknown",
    "Unknown",
    "Damage by reactive oxygen species",
    "Unknown",
    "Concurrent POLD1 mutations and defective DNA mismatch repair",
    "Defective DNA mismatch repair",
    "Aristolochic acid exposure",
    "Unknown",
    "Aflatoxin exposure",
    "Chemotherapy treatment",
    "Defective DNA mismatch repair",
    "Possible sequencing artefact",
    "Unknown",
    "Tobacco chewing",
    "Defective DNA base excision repair due to NTHL1 mutations",
    "Platinum chemotherapy treatment",
    "Azathioprine treatment",
    "Unknown",
    "Unknown",
    "Platinum chemotherapy treatment",
    "Defective DNA base excision repair due to MUTYH mutations",
    "Unknown",
    "Indirect effect of ultraviolet light",
    "Unknown",
    "Unknown",
    "Unknown",
    "Haloalkane exposure",
    "Possible sequencing artefact",
    "Defective DNA mismatch repair",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Possible sequencing artefact",
    "Activity of activation-induced cytidine deaminase (AID)",
    "Indirect effects of activation-induced cytidine deaminase (AID)",
    "Unknown chemotherapy treatment",
    "Thiopurine chemotherapy treatment",
    "Colibactin exposure (E.coli bacteria carrying pks pathogenicity island)",
    "Unknown",
    "Duocarmycin exposure"
  )
)

rownames(SBS_aetiology) <- colnames(SBS)

SBS <- apply(SBS, 2, function(x) x / sum(x)) %>% as.data.frame()

SBS_db <- list(
  db = SBS,
  aetiology = SBS_aetiology,
  date = "2020/09/29"
)

saveRDS(SBS_db, file = "inst/extdata/SBS_signatures.rds")

## DBS
DBS <- readr::read_csv("data-raw/SigProfiler Signatures/SigProfiler reference signatures/SigProfiler reference whole-genome signatures/sigProfiler_DBS_signatures.csv")
##
DBS_db <- list()
DBS <- DBS %>% column_to_rownames("Mutation Type")
## Ref link: https://cancer.sanger.ac.uk/cosmic/signatures/DBS
DBS_aetiology <- data.frame(
  aetiology = c(
    "Ultraviolet light exposure",
    "Tobacco smoking and other mutagens",
    "Polymerase epsilon exonuclease domain mutations",
    "Unknown",
    "Platinum chemotherapy treatment",
    "Unknown",
    "Defective DNA mismatch repair",
    "Unknown",
    "Unknown",
    "Defective DNA mismatch repair",
    "Unknown (possibly related to APOBEC mutagenesis)"
  )
)

rownames(DBS_aetiology) <- colnames(DBS)

DBS <- apply(DBS, 2, function(x) x / sum(x)) %>% as.data.frame()

DBS_db <- list(
  db = DBS,
  aetiology = DBS_aetiology,
  date = "2020/09/29"
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

cosmic_v3.1 <- readxl::read_excel("data-raw/COSMIC_Mutational_Signatures_v3.1.xlsx", sheet = 3)
cosmic_v3.1 <- cosmic_v3.1 %>% column_to_rownames("Type")
rownames(cosmic_v3.1)

cosmic_v3.1 <- cosmic_v3.1[rownames(ID), ]
get_sig_similarity(cosmic_v3.1, ID)

ID <- cosmic_v3.1
colSums(ID)

## Ref link: https://cancer.sanger.ac.uk/cosmic/signatures/ID
ID_aetiology <- data.frame(
  aetiology = c(
    "Slippage during DNA replication of the replicated DNA strand",
    "Slippage during DNA replication of the template DNA strand",
    "Tobacco smoking",
    "Unknown",
    "Unknown",
    "Defective homologous recombination DNA damage repair",
    "Defective DNA mismatch repair",
    "Repair of DNA double strand breaks by NHEJ mechanisms or mutations in topoisomerase TOP2A",
    "Unknown",
    "Unknown",
    "Unknown",
    "Unknown",
    "Ultraviolet light exposure",
    "Unknown",
    "Unknown",
    "Unknown",
    "Mutations in topoisomerase TOP2A",
    "Colibactin exposure (E.coli bacteria carrying pks pathogenicity island)"
  )
)

rownames(ID_aetiology) <- colnames(ID)
ID <- apply(ID, 2, function(x) x / sum(x)) %>% as.data.frame()

ID_db <- list(
  db = ID,
  aetiology = ID_aetiology,
  date = "2020/09/29"
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
    ),
    sep = ":"
  ))

TSB_db <- list()
TSB <- TSB %>%
  dplyr::select(-c("Strand", "Type", "Subtype")) %>%
  column_to_rownames("component")
## Ref link: https://cancer.sanger.ac.uk/cosmic/signatures/SBS
sbs_file <- "inst/extdata/SBS_signatures.rds"
sbs <- readRDS(sbs_file)
rbind(colnames(sbs$db[, 1:65]), colnames(TSB))

colnames(TSB) <- colnames(TSB) %>% str_remove("-E")

TSB_aetiology <- sbs$aetiology[rownames(sbs$aetiology) %in% colnames(TSB), , drop = FALSE]

TSB <- apply(TSB, 2, function(x) x / sum(x)) %>% as.data.frame()

TSB_db <- list(
  db = TSB,
  aetiology = TSB_aetiology,
  date = "2020/09/29"
)

saveRDS(TSB_db, file = "inst/extdata/TSB_signatures.rds")

## SBS for different genomes
sbs_list <- list()
sbs_list$hg19 <- readxl::read_excel("data-raw/SBS_signatures_genome_builds.xlsx") %>%
  dplyr::mutate(
    Type = paste0(
      substr(Subtype, 1, 1),
      "[",
      Type,
      "]",
      substr(Subtype, 3, 3)
    )
  ) %>%
  dplyr::select(-Subtype) %>%
  tibble::column_to_rownames("Type")

sbs_list$hg38 <- readxl::read_excel("data-raw/SBS_signatures_genome_builds.xlsx", sheet = 2) %>%
  dplyr::mutate(
    Type = paste0(
      substr(Subtype, 1, 1),
      "[",
      Type,
      "]",
      substr(Subtype, 3, 3)
    )
  ) %>%
  dplyr::select(-Subtype) %>%
  tibble::column_to_rownames("Type")
sbs_list$mm9 <- readxl::read_excel("data-raw/SBS_signatures_genome_builds.xlsx", sheet = 3) %>%
  dplyr::mutate(
    Type = paste0(
      substr(Subtype, 1, 1),
      "[",
      Type,
      "]",
      substr(Subtype, 3, 3)
    )
  ) %>%
  dplyr::select(-Subtype) %>%
  tibble::column_to_rownames("Type")
sbs_list$mm10 <- readxl::read_excel("data-raw/SBS_signatures_genome_builds.xlsx", sheet = 4) %>%
  dplyr::mutate(
    Type = paste0(
      substr(Subtype, 1, 1),
      "[",
      Type,
      "]",
      substr(Subtype, 3, 3)
    )
  ) %>%
  dplyr::select(-Subtype) %>%
  tibble::column_to_rownames("Type")

sbs_file <- "inst/extdata/SBS_signatures.rds"
sbs <- readRDS(sbs_file)
sbs_list$aetiology = sbs$aetiology[rownames(sbs$aetiology) %in% colnames(sbs_list$hg19), , drop = FALSE]
sbs_list$date = "2020/09/29"

saveRDS(sbs_list, file = "inst/extdata/SBS_signatures_list.rds")

## DBS for different genomes
dbs_list <- list()
dbs_list$hg19 <- readxl::read_excel("data-raw/DBS_signatures_genome_builds.xlsx") %>%
  tibble::column_to_rownames("Type")
dbs_list$hg38 <- readxl::read_excel("data-raw/DBS_signatures_genome_builds.xlsx", sheet = 2) %>%
  tibble::column_to_rownames("Type")
dbs_list$mm9 <- readxl::read_excel("data-raw/DBS_signatures_genome_builds.xlsx", sheet = 3) %>%
  tibble::column_to_rownames("Type")
dbs_list$mm10 <- readxl::read_excel("data-raw/DBS_signatures_genome_builds.xlsx", sheet = 4) %>%
  tibble::column_to_rownames("Type")

dbs_file <- "inst/extdata/DBS_signatures.rds"
dbs <- readRDS(dbs_file)
dbs_list$aetiology = dbs$aetiology
dbs_list$date = "2020/09/29"

saveRDS(dbs_list, file = "inst/extdata/DBS_signatures_list.rds")
