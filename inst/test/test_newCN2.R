load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))

# debug(sig_tally)
# cn_tally_W <- sig_tally(cn, method = "W")

cnlist <- sigminer:::get_cnlist(cn, add_index = TRUE)

features <- c("BP10MB", "CNCP", "OsCN")
for (feature in features) {
  if (feature == "BP10MB") {
    getPC_BP10MB(cnlist)
  }else if(feature == "CNCP"){
    getPC_CNCP_Max(cnlist)
  }else if(feature == "OsCN"){
    getPC_OsCN(cnlist)
  }
}

# abs_profiles <- cnlist
get_PC10MB <- function(abs_profiles) {
  y = purrr::map_df(abs_profiles, function(x){

    x_cp = data.table::copy(x)
    x_cp[, `:=`(
      region_start = start - 5e+06,
      region_end = start + 5e+06
    )][, region_start := ifelse(region_start <= 0, 1, region_start)]
    x_cp$start = NULL
    x_cp$end = NULL
    x_cp$segVal = NULL

    setkey(x_cp, "chromosome", "region_start", "region_end")
    x_overlap <- data.table::foverlaps(x, x_cp,
                                       by.x = c("chromosome", "start", "end"),
                                       type = "any")
    x %>%
      dplyr::as_tibble() %>%
      dplyr::full_join(
        x_overlap %>%
          dplyr::as_tibble() %>%
          dplyr::group_by(.data$Index) %>%
          dplyr::summarise(
            value = dplyr::n() - 1L
          ),
        by = c("Index")
      ) %>%
      dplyr::select(c("sample", "value", "Index"))
  }) %>%
    data.table::as.data.table()

  y[order(y$Index)]

}

getPC_CNCP_Left <- function(abs_profiles) {
  y = purrr::map_df(abs_profiles, function(x){

    x_cp = data.table::copy(x)

    x_cp <- x_cp %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(chromosome)
    x_cp <- x_cp %>% dplyr::mutate(value = abs(c(0, diff(segVal)))) %>%
      dplyr::select(c("sample", "value", "Index"))
  }) %>%
    data.table::as.data.table()
  y[order(y$Index)]

}
# c(0, diff(1:3))
# rev(c(0, diff(rev(1:3))))

getPC_CNCP_Right <- function(abs_profiles){
  y = purrr::map_df(abs_profiles, function(x){
    x_cp = data.table::copy(x)

    x_cp <- x_cp %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(chromosome)
    x_cp <- x_cp %>%
      dplyr::mutate(value = abs(rev(c(0, diff(rev(segVal)))))) %>%
      dplyr::select(c("sample", "value", "Index"))
  }) %>%
    data.table::as.data.table()
  y[order(y$Index)]
}

# getPC_CNCP_Max <- function(abs_profiles){
#   y = purrr::map_df(abs_profiles, function(x){
#     Right_df <- getPC_CNCP_Right(abs_profiles)
#     Left_df <- getPC_CNCP_Left(abs_profiles)
#     max_df <- Right_df
#
#   }) %>%
#     data.table::as.data.table()
#   y[order(y$Index)]
# }

getPC_OsCN <- function(){

}



# sigminer:::getChangepointCN(cnlist)

# purely count
## outputing columns "ID", "value"
getPC_BP10MB
getPC_CNCP_Left
getPC_CNCP_Right
getPC_CNCP_Max
# getPC_CNCP_type
getPC_OsCN






