## code to prepare `Nature-Cancer-Sigs` dataset goes here
# Data form Nik-Zainal lab work

library(readxl)
library(magrittr)

s1 <- read_excel("data-raw/NC-Ref-Sigs.xlsx")
s2 <- read_excel("data-raw/NC-Ref-Sigs.xlsx", sheet = 2)
s3 <- read_excel("data-raw/NC-Ref-Sigs.xlsx", sheet = 3)
s4 <- read_excel("data-raw/NC-Ref-Sigs.xlsx", sheet = 4)


# S1 and S2 ---------------------------------------------------------------

s1 <- s1 %>% tibble::column_to_rownames("Type") %>% as.matrix()
s1_names <- colnames(s1)
s1_map <- stringr::str_split(s1_names, pattern = " ", simplify = T) %>%
  as.data.frame() %>%
  dplyr::mutate(V2 = stringr::str_remove(V2, "\\("),
                V2 = stringr::str_remove(V2, "\\)"),
                V2 = stringr::str_remove(V2, "^[^_]+_"),
                V2 = stringr::str_remove(V2, "SoftTissue_"),
                V2 = stringr::str_remove(V2, "neck_"),
                V2 = stringr::str_replace_all(V2, "_", ","))

colnames(s1) <- s1_map$V1

SBS_Nik_lab_Organ <- list(
  db = s1,
  aetiology = s1_map %>% tibble::column_to_rownames("V1") %>% setNames("aetiology"),
  date = "2020/10/20"
)

saveRDS(SBS_Nik_lab_Organ, file = "inst/extdata/SBS_signatures_Nik_lab_Organ.rds")

s2 <- s2 %>% tibble::column_to_rownames("Type") %>% as.matrix()
s2_names <- colnames(s2)
s2_map <- stringr::str_split(s2_names, pattern = " ", simplify = T) %>%
  as.data.frame() %>%
  dplyr::mutate(V2 = stringr::str_remove(V2, "\\("),
                V2 = stringr::str_remove(V2, "\\)"),
                V2 = stringr::str_remove(V2, "^[^_]+_"),
                V2 = stringr::str_remove(V2, "SoftTissue_"),
                V2 = stringr::str_remove(V2, "neck_"),
                V2 = stringr::str_replace_all(V2, "_", ","))

colnames(s2) <- s2_map$V1

RS_Nik_lab_Organ <- list(
  db = s2,
  aetiology = s2_map %>% tibble::column_to_rownames("V1") %>% setNames("aetiology"),
  date = "2020/10/20"
)

saveRDS(RS_Nik_lab_Organ, file = "inst/extdata/RS_signatures_Nik_lab_Organ.rds")

# S3 and S4 ---------------------------------------------------------------

s3 <- s3 %>% tibble::column_to_rownames("Type") %>% as.matrix()
colnames(s3) <- gsub("\\.", "_", colnames(s3))

SBS_Nik_lab <- list(
  db = s3,
  aetiology = data.frame(
    V1 = colnames(s3),
    V2 = "See COSMIC signatures with same ID or https://doi.org/10.1038/s43018-020-0027-5"
  ) %>% tibble::column_to_rownames("V1") %>% setNames("aetiology"),
  date = "2020/10/20"
)

saveRDS(SBS_Nik_lab, file = "inst/extdata/SBS_signatures_Nik_lab.rds")

s4 <- s4 %>% tibble::column_to_rownames("Type") %>% as.matrix()
colnames(s4) <- gsub("\\.", "_", colnames(s4))

RS_Nik_lab <- list(
  db = s4,
  aetiology = data.frame(
    V1 = colnames(s4),
    V2 = "See https://doi.org/10.1038/s43018-020-0027-5"
  ) %>% tibble::column_to_rownames("V1") %>% setNames("aetiology"),
  date = "2020/10/20"
)

saveRDS(RS_Nik_lab, file = "inst/extdata/RS_signatures_Nik_lab.rds")
