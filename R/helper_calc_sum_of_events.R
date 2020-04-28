calculateSumOfEvents <-
  function(feature,
           component,
           name,
           type) {
    params <- flexmix::parameters(component)
    if (!is.null(nrow(params))) {
      comp_orders <- order(params[1, ])
    } else {
      comp_orders <- order(params)
    }

    if (type == "count") {
      df <- feature
      send_success("Calculating the sum of event counts for ", name, ".")
      df$cluster <- flexmix::clusters(component, data.frame(dat = as.numeric(feature[, 2]))) %>%
        as.character()
      df$cluster <- factor(df$cluster, levels = seq_along(comp_orders) %>% as.character())
      df$count <- 1
      df_sum <- df %>%
        dplyr::select(-.data$value) %>%
        dplyr::group_by(.data$ID, .data$cluster) %>%
        dplyr::summarise(count = sum(.data$count)) %>%
        tidyr::spread(key = "cluster", value = "count", fill = 0, drop = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(c("ID", as.character(comp_orders)))
    } else {
      send_info("Calculating the sum of posterior probabilities for ", name, ".")
      curr <- flexmix::posterior(component, data.frame(dat = as.numeric(feature[, 2])))
      df <- cbind(feature, curr)
      df_sum <- df %>%
        dplyr::select(-.data$value) %>%
        dplyr::group_by(.data$ID) %>%
        dplyr::summarise_if(is.numeric, sum) %>%
        dplyr::select(c("ID", as.character(comp_orders)))
    }

    mat <- df_sum %>%
      dplyr::select(-.data$ID) %>%
      as.matrix()

    rownames(mat) <- df_sum$ID
    colnames(mat) <- paste0(name, 1:ncol(mat))

    mat
  }
