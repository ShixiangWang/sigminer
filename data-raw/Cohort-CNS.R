# USARC CNS ----------------------------------------------------------------

usarc <- read_excel("data-raw/USARC.CNS.xlsx", sheet = 1, skip = 1)
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
  tidyr::unite(col = "class", c("cluster", "type", "size"), sep = "_") %>%
  dplyr::mutate(class = gsub("_NA", "", class)) %>%
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
