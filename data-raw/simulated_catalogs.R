## code to prepare `simulated_catalogs` dataset goes here
# data from https://doi.org/10.1038/s43018-020-0027-5
# 5 simulated catalog datasets are used by the paper but only 4 are available
data1 <- readxl::read_excel("data-raw/NC-Supp-data.xlsx", sheet = 1)
data2 <- readxl::read_excel("data-raw/NC-Supp-data.xlsx", sheet = 2)
data3 <- readxl::read_excel("data-raw/NC-Supp-data.xlsx", sheet = 3)
data4 <- readxl::read_excel("data-raw/NC-Supp-data.xlsx", sheet = 4)
# data5 <- readxl::read_excel("data-raw/NC-Supp-data.xlsx", sheet = 5)

simulated_catalogs <- list(
  data1, data2, data3, data4
) %>%
  purrr::map(.f = function(x) {
    x %>%
      tibble::column_to_rownames("...1") %>%
      as.matrix()
  }) %>%
  setNames(c("set1", "set2", "set3", "set4"))

usethis::use_data(simulated_catalogs, overwrite = TRUE)
