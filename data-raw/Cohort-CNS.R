# USARC CNS ----------------------------------------------------------------

usarc <- read_excel("data-raw/USARC.CNS.xlsx", sheet = 1, skip = 1)
colnames(usarc)[1] <- "class"

CN40_map <- usarc %>%
  dplyr::select("class") %>%
  tidyr::separate("class", c("type", "type2", "size"), sep = ":", remove = FALSE) %>%
  tidyr::separate("size", c("size_min", "size_max"), sep = ",") %>%
  dplyr::mutate(
    size_min = as.numeric(sub("\\(", "", size_min)),
    size_max = as.numeric(sub("\\]", "", size_max))
  ) %>%
  dplyr::mutate(
    type = dplyr::case_when(
      type == "amp" ~ "5+",
      type == "dup" ~ "3-4",
      type == "neut" ~ "2",
      type == "del" ~ "0-1"
    ),
    size = dplyr::case_when(
      size_min < 0 & size_max < 1 ~ paste0("0-", 1000*size_max, "Kb"),
      #size_min < 0 & size_max >= 1 ~ paste0("0-", size_max, "Mb")
      is.infinite(size_max) ~ paste0(">", size_min, "Mb"),
      size_max < 1 ~ paste0(1000*size_min, "Kb-", 1000*size_max, "Kb"),
      size_min < 1 & size_max == 1 ~ paste0(1000*size_min, "Kb-", size_max, "Mb"),
      size_min >= 1 ~ paste0(size_min, "Mb-", size_max, "Mb")
    )
  ) %>%
  tidyr::unite("label", c("type", "type2", "size"), sep = ":") %>%
  dplyr::mutate(
    cn = as.integer(substr(label, 1, 1))
  ) %>%
  dplyr::arrange(cn) %>%
  dplyr::select(class, label) %>%
  dplyr::rename(
    label1 = class,
    label2 = label
  )

data.table::fwrite(CN40_map, file = "inst/extdata/CN40-Map.txt", sep = "\t")


map <- CN40_map$label2
names(map) <- CN40_map$label1

usarc2 <- usarc %>%
  tibble::column_to_rownames("class")
rownames(usarc2) <- map[rownames(usarc2)]

colnames(usarc2) <- paste0("USARC_", colnames(usarc2))

# check sum
colSums(usarc2)

CNS_USARC <- list(
  db = as.matrix(usarc2),
  aetiology = data.frame(
    V1 = colnames(usarc2),
    V2 = "See https://doi.org/10.1016/j.ccell.2019.02.002"
  ) %>% tibble::column_to_rownames("V1") %>% setNames("aetiology"),
  date = "2021/01/15"
)

saveRDS(CNS_USARC, file = "inst/extdata/CNS_signatures_USARC.rds")
