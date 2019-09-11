calculateSumOfPosteriors <-
  function(CN_feature,
             components,
             name,
             rowIter = 1000,
             cores = 1) {
    if (cores > 1) {

      # attachNamespace("foreach")
      # attachNamespace("doParallel")

      len <- dim(CN_feature)[1]
      iters <- len %/% rowIter
      if (len %% rowIter == 0) {
        iters <- iters - 1
      }

      ## Qs: Should put the parallel out of this fun?
      doFuture::registerDoFuture()
      oplan <- future::plan()
      future::plan("multisession", workers = cores)
      on.exit(future::plan(oplan), add = TRUE)

      curr_posterior <- foreach::foreach(i = 0:iters, .combine = rbind) %dopar% {
        start <- i * rowIter + 1
        if (i != iters) {
          end <- (i + 1) * rowIter
        } else {
          end <- len
        }
        flexmix::posterior(components, data.frame(dat = as.numeric(CN_feature[start:end, 2])))
      }
    } else {
      curr_posterior <-
        flexmix::posterior(components, data.frame(dat = as.numeric(CN_feature[, 2])))
    }

    mat <- cbind(CN_feature, curr_posterior)
    posterior_sum <- c()

    ## foreach and parallelising doesn't make the following code faster.
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
