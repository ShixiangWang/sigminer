#' Extract Signatures with SigProfiler
#'
#' This function provides an interface to software SigProfiler.
#' More please see <https://github.com/AlexandrovLab/SigProfilerExtractor>.
#' Typically, a reference genome is not required because the input is a matrix (my understanding).
#' **If you are using refitting result by SigProfiler, please make sure you have input the matrix same order as examples at <https://github.com/AlexandrovLab/SigProfilerMatrixGenerator/tree/master/SigProfilerMatrixGenerator/references/matrix/BRCA_example>**. If not, use `sigprofiler_reorder()` firstly.
#'
#' @inheritParams sig_extract
#' @name sigprofiler
#' @rdname sigprofiler
#' @param output output directory.
#' @param output_matrix_only if `TRUE`, only generate matrix file for SigProfiler
#' so user can call SigProfiler with the input by himself.
#' @param range signature number range, i.e. `2:5`.
#' @param nrun the number of iteration to be performed to extract each signature number.
#' @param refit if `TRUE`, then refit the denovo signatures with nnls. Same
#' meaning as `optimize` option in [sig_extract] or [sig_auto_extract].
#' @param refit_plot if `TRUE`, SigProfiler will make
#' denovo to COSMIC sigantures decompostion plots. However, this may fail due
#' to some matrix cannot be identified by SigProfiler plot program.
#' @param is_exome if `TRUE`, the exomes will be extracted.
#' @param init_method the initialization algorithm for W and H matrix of NMF.
#' Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar', 'alexandrov-lab-custom'
#' and 'nndsvd_min'.
#' @param cores number of cores used for computation.
#' @param genome_build I think this option is useless when input is `matrix`, keep it
#' in case it is useful.
#' @param use_conda if `TRUE`, create an independent conda environment to run SigProfiler.
#' @param py_path path to Python executable file, e.g. '/Users/wsx/anaconda3/bin/python'.
#' @param sigprofiler_version version of `SigProfilerExtractor`. If this
#' package is not installed, the specified package will be installed.
#' If this package is installed, this option is useless.
#'
#' @return For `sigprofiler_extract()`, returns nothing. See `output` directory.
#' @export
#'
#' @examples
#' if (FALSE) {
#'   load(system.file("extdata", "toy_copynumber_tally_W.RData",
#'     package = "sigminer", mustWork = TRUE
#'   ))
#'
#'   reticulate::conda_list()
#'
#'   sigprofiler_extract(cn_tally_W$nmf_matrix, "~/test/test_sigminer",
#'     use_conda = TRUE
#'   )
#'
#'   sigprofiler_extract(cn_tally_W$nmf_matrix, "~/test/test_sigminer",
#'     use_conda = FALSE, py_path = "/Users/wsx/anaconda3/bin/python"
#'   )
#' }
sigprofiler_extract <- function(nmf_matrix, output,
                                output_matrix_only = FALSE,
                                range = 2:5, nrun = 10L,
                                refit = FALSE,
                                refit_plot = FALSE,
                                is_exome = FALSE,
                                init_method = c(
                                  "random", "nndsvd_min",
                                  "nndsvd", "nndsvda", "nndsvdar"
                                ),
                                cores = -1L,
                                genome_build = c("hg19", "hg38", "T2T", "mm10", "mm9", "ce11"),
                                use_conda = FALSE,
                                py_path = NULL,
                                sigprofiler_version = "1.1.3") {
  output <- path.expand(output)
  genome_build <- match.arg(genome_build)

  genome_build <- switch(genome_build,
    hg19 = "GRCh37",
    hg38 = "GRCh38",
    genome_build
  )

  if (!dir.exists(output)) {
    dir.create(output, recursive = TRUE)
  }

  env_install(use_conda, py_path, pkg = "SigProfilerExtractor", pkg_version = sigprofiler_version)

  init_method <- match.arg(init_method)
  cores <- as.integer(cores)
  nrun <- as.integer(nrun)

  sig_ranges <- range(range)
  sys <- reticulate::import("sys")
  sigpro <- reticulate::import("SigProfilerExtractor.sigpro")

  # print(sys$path)

  nmf_matrix <- t(nmf_matrix)

  ii <- colSums(nmf_matrix) < 0.01
  if (any(ii)) {
    message(
      "The follow samples dropped due to null catalogue:\n\t",
      paste0(colnames(nmf_matrix)[ii], collapse = ", ")
    )
    nmf_matrix <- nmf_matrix[, !ii, drop = FALSE]
  }

  in_df <- nmf_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Mutation Types")

  if (output_matrix_only) {
    if (!dir.exists(output)) dir.create(output, recursive = TRUE)
    data.table::fwrite(data.table::as.data.table(in_df),
                       file = file.path(output, "sigprofiler_input.txt"),
                       sep = "\t"
    )
    message("SigProfiler input file saved to ", file.path(output, "sigprofiler_input.txt"))
    return(NULL)
  }

  tmp_dir <- basename(tempfile(pattern = "dir"))
  tmp_dir_full <- file.path(tempdir(), tmp_dir)
  if (dir.exists(tmp_dir_full)) {
    unlink(tmp_dir_full)
  }
  dir.create(tmp_dir_full, recursive = TRUE)
  tmp_file <- file.path(tmp_dir_full, "sigprofiler_input.txt")
  data.table::fwrite(data.table::as.data.table(in_df),
    file = tmp_file,
    sep = "\t"
  )
  on.exit(unlink(tmp_file))

  ## The newst version
  ## https://github.com/AlexandrovLab/SigProfilerExtractor
  tryCatch(
    {
      sigpro$sigProfilerExtractor(
        "matrix",
        output,
        tmp_file,
        reference_genome = genome_build,
        opportunity_genome = genome_build,
        minimum_signatures = sig_ranges[1],
        maximum_signatures = sig_ranges[2],
        nmf_replicates = nrun,
        exome = is_exome,
        nmf_init = init_method,
        refit_denovo_signatures = refit,
        make_decomposition_plots = refit_plot,
        cpu = cores
      )
      sys$stdout$flush()
    },
    error = function(e) {
      cat(e$message)
      message("Error detected, try recalling with python command.")
      if (is.null(py_path)) {
        message("Python path not set, try calling with 'python'.")
        py_path <- "python"
      }

      cmd_import <- "from SigProfilerExtractor import sigpro as sig"
      cmd_call <- paste(
        quote_opt("matrix"),
        quote_opt(output),
        quote_opt(tmp_file),
        quote_opt(genome_build, opt = "reference_genome"),
        quote_opt(genome_build, opt = "opportunity_genome"),
        quote_opt(sig_ranges[1], opt = "minimum_signatures", rm_quote = TRUE),
        quote_opt(sig_ranges[2], opt = "maximum_signatures", rm_quote = TRUE),
        quote_opt(nrun, opt = "nmf_replicates", rm_quote = TRUE),
        quote_opt(ifelse(is_exome, "True", "False"), opt = "exome", rm_quote = TRUE),
        quote_opt(init_method, opt = "nmf_init"),
        quote_opt(ifelse(refit, "True", "False"), opt = "refit_denovo_signatures", rm_quote = TRUE),
        quote_opt(ifelse(refit_plot, "True", "False"), opt = "make_decomposition_plots", rm_quote = TRUE),
        quote_opt(cores, opt = "cpu", rm_quote = TRUE),
        sep = ","
      )

      cmds <- paste(
        py_path, "-c",
        paste0(
          "\"",
          cmd_import,
          ";",
          "sig.sigProfilerExtractor(",
          cmd_call,
          ")\""
        )
      )

      message("===")
      message(cmds)
      message("===")
      system(cmds)
    }
  )
}

quote_opt <- function(value, opt = NULL, rm_quote = FALSE) {
  if (is.null(opt)) {
    if (rm_quote) {
      value
    } else {
      paste0("'", value, "'")
    }
  } else {
    if (rm_quote) {
      paste(opt,
        value,
        sep = "="
      )
    } else {
      paste(opt,
        paste0("'", value, "'"),
        sep = "="
      )
    }
  }
}



# Import SigProfiler Results ----------------------------------------------

#' Import SigProfiler Results into R
#'
#' @inheritParams sigprofiler
#' @param order_by_expo if `TRUE`, order the import signatures by their exposures, e.g. the signature
#' contributed the most exposure in all samples will be named as `Sig1`.
#' @param type one of 'suggest' (for suggested solution), 'refit' (for refit solution) or 'all' (for all solutions).
#'
#' @return For `sigprofiler_import()`, a `list` containing `Signature` object.
#' @export
#' @rdname sigprofiler
sigprofiler_import <- function(output,
                               order_by_expo = FALSE,
                               type = c("suggest", "refit", "all")) {
  stopifnot(dir.exists(output))
  type <- match.arg(type)

  result_dir <- list.dirs(output, full.names = TRUE, recursive = FALSE)

  if (length(result_dir) != 1) {
    if (length(result_dir) == 0) {
      stop("Cannot find any result directory in ", output, "!")
    } else {
      stop("Should not have more than 1 directory be found, please check ", output, "!")
    }
  }

  all_solution_stats <- data.table::fread(file.path(result_dir, "All_solutions_stat.csv"))
  colnames(all_solution_stats)[1] <- "SignatureNumber"

  message("NOTE: signature(A,B,C)... will be renamed to Sig(1,2,3)...")

  if (type %in% c("suggest", "refit")) {
    solution_path <- file.path(result_dir, "Suggested_Solution")
    if (type == "suggest") {
      solution_path <- list.files(solution_path, pattern = "Novo", full.names = TRUE)
      message("Reading suggested solution...")
    } else {
      solution_path <- list.files(solution_path, pattern = "Decomposed", full.names = TRUE)
      message("Reading suggested solution...")
    }

    if (length(solution_path) != 1) {
      stop("No solution path or more than 1 solution path found, please check!")
    }

    solution <- read_sigprofiler_solution(solution_path, order_by_expo = order_by_expo)

    message("Done.")
    return(list(
      solution = solution,
      all_stats = all_solution_stats
    ))
  } else {
    solution_path <- file.path(result_dir, "All_Solutions")
    message("Reading all solutions...")
    solutions_path <- list.dirs(solution_path, full.names = TRUE, recursive = FALSE)
    solutions <- purrr::map(solutions_path, read_sigprofiler_solution, order_by_expo = order_by_expo)
    names(solutions) <- paste0("S", sub("[^_]+_(.+)_[^_]+", "\\1", basename(solutions_path)))

    message("Done.")
    return(list(
      solution_list = solutions,
      all_stats = all_solution_stats
    ))
  }
}

read_sigprofiler_solution <- function(x, order_by_expo = FALSE) {
  expo_path <- list.files(x, pattern = "Activities.*.txt", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(expo_path) > 1) {
    expo_path <- expo_path[!grepl("error", expo_path, ignore.case = TRUE)]
  }
  sigs_path <- list.files(x, pattern = "Signatures.txt", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  stat_samp_path <- list.files(x, pattern = "Samples_Stats.*.txt", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  stat_sigs_path <- list.files(x, pattern = "Signatures_Stats.txt", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(stat_sigs_path) < 1) {
    message("No Signatures_Stats.txt found, try finding *map_to_COSMIC* file.")
    stat_sigs_path <- list.files(x, pattern = "map_to_COSMIC", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    refit <- TRUE
  } else {
    refit <- FALSE
  }

  expo <- data.table::fread(expo_path, header = TRUE)
  sigs <- data.table::fread(sigs_path, header = TRUE)
  stat_samp <- data.table::fread(stat_samp_path, header = TRUE)
  stat_sigs <- data.table::fread(stat_sigs_path, header = TRUE)

  K <- ncol(expo) - 1L

  colnames(expo) <- c("sample", paste0("Sig", seq_len(K)))
  colnames(sigs) <- c("component", paste0("Sig", seq_len(K)))
  if (refit) {
    colnames(stat_sigs)[2] <- "Signatures"
  } else {
    colnames(stat_sigs)[1] <- "Signatures"
  }
  stat_sigs$Signatures <- colnames(expo)[-1]
  colnames(stat_samp)[1] <- "Samples"

  Signature.norm <- sigs %>%
    dplyr::as_tibble() %>%
    tibble::column_to_rownames("component") %>%
    as.matrix()

  Exposure <- expo %>%
    dplyr::as_tibble() %>%
    tibble::column_to_rownames("sample") %>%
    as.matrix() %>%
    t()

  Exposure.norm <- apply(Exposure, 2, function(x) x / sum(x, na.rm = TRUE))
  # When only one signature
  if (!is.matrix(Exposure.norm)) {
    Exposure.norm <- matrix(Exposure.norm, nrow = 1, dimnames = list(
      rownames(Exposure),
      colnames(Exposure)
    ))
  }

  Signature <- Signature.norm
  for (j in seq_len(K)) {
    Signature[, j] <- Signature[, j] * rowSums(Exposure)[j]
  }

  if (order_by_expo) {
    sig_orders <- sort(rowSums(Exposure), decreasing = TRUE) %>% names()
    if (!identical(rownames(Exposure), sig_orders)) {
      Exposure <- Exposure[sig_orders, , drop = FALSE]
      Exposure.norm <- Exposure.norm[sig_orders, , drop = FALSE]
      Signature <- Signature[, sig_orders, drop = FALSE]
      Signature.norm <- Signature.norm[, sig_orders, drop = FALSE]
      rownames(Exposure) <- rownames(Exposure.norm) <- colnames(Signature) <- colnames(Signature.norm) <- paste0("Sig", seq_along(sig_orders))
      stat_sigs$Signatures <- sig_orders
    }
  }

  res <- list(
    Signature = Signature,
    Signature.norm = Signature.norm,
    Exposure = Exposure,
    Exposure.norm = Exposure.norm,
    K = K,
    Stats = list(
      samples = stat_samp,
      signatures = stat_sigs
    )
  )

  class(res) <- "Signature"
  attr(res, "call_method") <- "SigProfiler"

  res
}


# Reorder matrix orders ---------------------------------------------------

#' Reorder Mutational Types Based on SigProfiler Reference Example File
#'
#' @inheritParams sig_extract
#' @param type mutational signature type.
#'
#' @return A NMF matrix for input of `sigprofiler_extract()`.
#' @export
#' @rdname sigprofiler
#'
#' @examples
#'
#' data("simulated_catalogs")
#' sigprofiler_reorder(t(simulated_catalogs$set1))
sigprofiler_reorder = function(nmf_matrix,
                               type = c("SBS96", "SBS6", "SBS12", "SBS192", "SBS1536", "SBS3072",
                                        "DBS78", "DBS312", "DBS1248", "DBS4992")
                               ) {
  type = match.arg(type)
  message("Downloading reference file...")
  ref = data.table::fread(sprintf("https://raw.githubusercontent.com/AlexandrovLab/SigProfilerMatrixGenerator/master/SigProfilerMatrixGenerator/references/matrix/BRCA_example/BRCA_example.%s.all", type))
  ref = tibble::column_to_rownames(ref, "MutationType")
  ref = t(as.matrix(ref))

  if (identical(colnames(nmf_matrix), colnames(ref))) {
    message("Order of your mutational types are correct, no need to reorder.")
  } else {
    message("Reordering...")
    nmf_matrix = nmf_matrix[, colnames(ref)]
    message("Done")
  }

  return(nmf_matrix)
}
