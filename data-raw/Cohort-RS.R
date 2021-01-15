library(readxl)

# BRCA RS -----------------------------------------------------------------
# Ref: https://www.nature.com/articles/nature17676
brca560 <- read_excel("data-raw/BRCA-RS.xlsx", sheet = 2)
brca560_2 <- brca560 %>%
  dplyr::mutate(
    class = paste(Type, Size, sep = ":"),
    class = gsub("k", "K", class)
  ) %>%
  dplyr::select(-Type, -Size) %>%
  tibble::column_to_rownames("class")

ref <- get_sig_db("RS_Nik_lab")

rownames(brca560_2)
rownames(ref$db)
which(!rownames(brca560_2) %in% rownames(ref$db))
nrow(brca560_2)
nrow(ref$db)

colnames(brca560_2) <- gsub("Rearrangement Signature ", "BRCA_RS", colnames(brca560_2))

# check sum
colSums(brca560_2)
brca560_2 <- apply(brca560_2, 2, function(x) x / sum(x))
colSums(brca560_2)

RS_brca560 <- list(
  db = brca560_2,
  aetiology = data.frame(
    V1 = colnames(brca560_2),
    V2 = "See https://doi.org/10.1038/nature17676"
  ) %>% tibble::column_to_rownames("V1") %>% setNames("aetiology"),
  date = "2021/01/15"
)

saveRDS(RS_brca560, file = "inst/extdata/RS_signatures_BRCA560.rds")

# USARC RS ----------------------------------------------------------------

usarc <- read_excel("data-raw/USARC.RS.xlsx", sheet = 1, skip = 1)
colnames(usarc)[1] <- "class"

usarc2 <- usarc %>%
  tidyr::separate(col = "class", into = c("type", "size", "cluster"), sep = ":") %>%
  dplyr::mutate(
    type = dplyr::case_when(
      type == "deletion" ~ "del",
      type == "inversion" ~ "inv",
      type == "tandem-duplication" ~ "tds",
      type == "translocation" ~ "trans"
    ),
    cluster = ifelse(cluster == "clustered", "clustered", "non-clustered"),
    size = sub("k", "K", size)
  ) %>%
  tidyr::unite(col = "class", c("cluster", "type", "size"), sep = ":") %>%
  dplyr::mutate(class = gsub(":NA", "", class)) %>%
  tibble::column_to_rownames("class")

rownames(usarc2)
rownames(ref$db)
which(!rownames(usarc2) %in% rownames(ref$db))
nrow(usarc2)
nrow(ref$db)

colnames(usarc2) <- gsub("\\.", "_", colnames(usarc2))

# check sum
colSums(usarc2)

RS_USARC <- list(
  db = as.matrix(usarc2),
  aetiology = data.frame(
    V1 = colnames(usarc2),
    V2 = "See https://doi.org/10.1016/j.ccell.2019.02.002"
  ) %>% tibble::column_to_rownames("V1") %>% setNames("aetiology"),
  date = "2021/01/15"
)

saveRDS(RS_USARC, file = "inst/extdata/RS_signatures_USARC.rds")
