#' Extract Signatures with SigProfiler
#'
#' This function provides an interface to software SigProfiler.
#' More please see <https://github.com/AlexandrovLab/SigProfilerExtractor>.
#' Typically, a reference genome is not required because the input is a matrix (my understanding).
#'
#' @inheritParams sig_extract
#' @name sigprofiler
#' @rdname sigprofiler
#' @param output output directory.
#' @param range signature number range, i.e. `2:5`.
#' @param nrun the number of iteration to be performed to extract each signature number.
#' @param is_exome if `TRUE`, the exomes will be extracted.
#' @param init_method the initialization algorithm for W and H matrix of NMF.
#' Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar' and 'alexandrov-lab-custom'.
#' @param cores number of cores used for computation.
#' @param genome_build I think this option is useless when input is `matrix`, keep it
#' in case it is useful.
#' @param use_conda if `TRUE`, create an independent conda environment to run SigProfiler.
#' @param py_path path to Python executable file, e.g. '/Users/wsx/anaconda3/bin/python'.
#' @param sigprofiler_version version of `SigProfilerExtractor`.
#'
#' @return For `sigprofiler_extract()`, returns nothing. See `output` directory.
#' @export
#'
#' @examples
#' if (FALSE) {
#'   load(system.file("extdata", "toy_copynumber_tally_M.RData",
#'     package = "sigminer", mustWork = TRUE
#'   ))
#'
#'   reticulate::conda_list()
#'
#'   sigprofiler_extract(cn_tally_M$nmf_matrix, "~/test/test_sigminer",
#'     use_conda = TRUE
#'   )
#'
#'   sigprofiler_extract(cn_tally_M$nmf_matrix, "~/test/test_sigminer",
#'     use_conda = FALSE, py_path = "/Users/wsx/anaconda3/bin/python"
#'   )
#' }
sigprofiler_extract <- function(nmf_matrix, output, range = 2:5, nrun = 10L,
                                is_exome = TRUE,
                                init_method = c("random", "alexandrov-lab-custom", "nndsvd", "nndsvda", "nndsvdar"),
                                cores = -1L,
                                genome_build = c("hg19", "hg38", "mm10"),
                                use_conda = FALSE,
                                py_path = NULL,
                                sigprofiler_version = "1.0.17") {
  output <- path.expand(output)
  genome_build <- match.arg(genome_build)

  genome_build <- switch(genome_build,
    hg19 = "GRCh37",
    hg38 = "GRCh38",
    mm10 = "mm10"
  )

  if (!dir.exists(output)) {
    dir.create(output, recursive = TRUE)
  }

  if (!requireNamespace("reticulate")) {
    stop("Package 'reticulate' is required, please install it firstly!")
  }

  if (use_conda) {
    tryCatch(reticulate::conda_binary(),
      error = function(e) {
        message("Cannot find conda binary, installing miniconda...")
        reticulate::install_miniconda()
      }
    )

    ## Prepare conda environment and packages
    tryCatch(
      reticulate::use_condaenv("sigminer_sigprofiler", required = TRUE),
      error = function(e) {
        message("Conda environment not detected, creat it and install required packages.")
        message("======")
        reticulate::conda_create("sigminer_sigprofiler")
        message("Installing packages, be patient...")
        message("======")
        reticulate::conda_install("sigminer_sigprofiler",
          packages = paste0("SigProfilerExtractor==", sigprofiler_version),
          pip = TRUE
        )
        reticulate::use_condaenv("sigminer_sigprofiler", required = TRUE)
      }
    )

    message("Python and conda environment configuration.")
    message("====================")
    print(reticulate::py_config())
  } else {
    if (is.null(py_path)) {
      if (!reticulate::py_available(initialize = TRUE)) {
        stop("Python is not available!")
      }
      config <- reticulate::py_config()
      print(config)
      reticulate::use_python(config$python)
    } else {
      reticulate::use_python(py_path, required = TRUE)
    }
  }

  if (!reticulate::py_module_available("SigProfilerExtractor")) {
    message("Python module 'SigProfilerExtractor' not found, try installing it...")
    reticulate::py_install(paste0("SigProfilerExtractor==", sigprofiler_version), pip = TRUE)
  }

  init_method <- match.arg(init_method)
  cores <- as.integer(cores)
  nrun <- as.integer(nrun)

  sig_ranges <- range(range)
  sys <- reticulate::import("sys")
  sigpro <- reticulate::import("SigProfilerExtractor.sigpro")

  # print(sys$path)

  in_df <- nmf_matrix %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("MutationType")

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
        reference_genome = genome_build, opportunity_genome = genome_build,
        minimum_signatures = sig_ranges[1],
        maximum_signatures = sig_ranges[2],
        nmf_replicates = nrun,
        exome = is_exome,
        nmf_init = init_method,
        refit_denovo_signatures = FALSE,
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
        quote_opt("False", opt = "refit_denovo_signatures", rm_quote = TRUE),
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
#' @param type one of 'suggest' (for suggested solution) or 'all' (for all solutions).
#'
#' @return For `sigprofiler_import()`, a `list` containing `Signature` object.
#' @export
#' @rdname sigprofiler
sigprofiler_import <- function(output, type = c("suggest", "all")) {
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

  if (type == "suggest") {
    solution_path <- file.path(result_dir, "Suggested_Solution")
    solution_path <- list.files(solution_path, pattern = "Novo", full.names = TRUE)
    message("Reading suggested solution...")

    solution <- read_sigprofiler_solution(solution_path)

    message("Done.")
    return(list(
      solution = solution,
      all_stats = all_solution_stats
    ))
  } else {
    solution_path <- file.path(result_dir, "All_Solutions")
    message("Reading all solutions...")
    solutions_path <- list.dirs(solution_path, full.names = TRUE, recursive = FALSE)
    solutions <- purrr::map(solutions_path, read_sigprofiler_solution)
    names(solutions) <- paste0("S", sub("[^_]+_(.+)_[^_]+", "\\1", basename(solutions_path)))

    message("Done.")
    return(list(
      solution_list = solutions,
      all_stats = all_solution_stats
    ))
  }
}

read_sigprofiler_solution <- function(x) {
  expo_path <- list.files(x, pattern = "Activities.*.txt", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(expo_path) > 1) {
    expo_path <- expo_path[!grepl("error", expo_path, ignore.case = TRUE)]
  }
  sigs_path <- list.files(x, pattern = "Signatures.txt", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  stat_samp_path <- list.files(x, pattern = "Samples_Stats.*.txt", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  stat_sigs_path <- list.files(x, pattern = "Signatures_Stats.txt", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

  expo <- data.table::fread(expo_path)
  sigs <- data.table::fread(sigs_path)
  stat_samp <- data.table::fread(stat_samp_path)
  stat_sigs <- data.table::fread(stat_sigs_path)

  K <- ncol(expo) - 1L

  colnames(expo) <- c("sample", paste0("Sig", seq_len(K)))
  colnames(sigs) <- c("component", paste0("Sig", seq_len(K)))
  colnames(stat_sigs)[1] <- "Signatures"
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
