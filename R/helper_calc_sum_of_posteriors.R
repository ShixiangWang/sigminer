calculateSumOfPosteriors <-
  function(CN_feature,
             components,
             name) {

    curr_posterior <-
      flexmix::posterior(components, data.frame(dat = as.numeric(CN_feature[, 2])))


    mat <- cbind(CN_feature, curr_posterior)
    posterior_sum <- c()

    ## TODO: I think I can speed this up in a tidy way
    for (i in unique(mat$ID))
    {
      posterior_sum <-
        rbind(posterior_sum, colSums(mat[mat$ID == i, c(-1, -2), drop = FALSE]))
    }
    params <- flexmix::parameters(components)
    if (!is.null(nrow(params))) {
      posterior_sum <- posterior_sum[, order(params[1, ]), drop = FALSE]
    }
    else {
      posterior_sum <- posterior_sum[, order(params), drop = FALSE]
    }
    colnames(posterior_sum) <-
      paste0(name, 1:ncol(posterior_sum))
    rownames(posterior_sum) <- rownames(unique(mat$ID))
    posterior_sum
  }
