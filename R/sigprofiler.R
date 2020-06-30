#' Extract Signatures with SigProfiler
#'
#' This function provides an interface to software SigProfiler.
#' More please see <https://github.com/AlexandrovLab/SigProfilerExtractor>.
#' Typically, a reference genome is not required because the input is a matrix (my understanding).
#'
#' @inheritParams sig_extract
#' @rdname sigprofiler
#' @param output output directory.
#' @param range signature number range, i.e. `2:5`.
#' @param nrun the number of iteration to be performed to extract each number signature.
#' @param is_exome if `TRUE`, the exomes will be extracted.
#' @param init_method the initialization algorithm for W and H matrix of NMF.
#' Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar' and 'alexandrov-lab-custom'.
#' @param cores number of cores used for computation.
#' @param genome_build I think this option is useless when input is `matrix`, keep it
#' in case it is useful.
#' @param use_conda if `TRUE`, create an independent conda environment to run SigProfiler.
#' @param py_path path to Python executable file, e.g. '/Users/wsx/anaconda3/bin/python'.
#'
#' @return Nothing. See `output` directory.
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
                                py_path = NULL) {
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
             })

    ## Prepare conda environment and packages
    tryCatch(
      reticulate::use_condaenv("sigminer_sigprofiler", required = TRUE),
      error = function(e) {
        message("Conda environment not detected, creat it and install required packages.")
        reticulate::conda_create("sigminer_sigprofiler")
        message("Be patient...")
        reticulate::conda_install("sigminer_sigprofiler",
          packages = "SigProfilerExtractor",
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

    if (!reticulate::py_module_available("SigProfilerExtractor")) {
      message("Python module 'SigProfilerExtractor' not found, try installing it...")
      reticulate::py_install("SigProfilerExtractor", pip = TRUE)
    }
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
  # sigProfilerExtractor(input_type, out_put, input_data, reference_genome="GRCh37",
  #                      opportunity_genome = "GRCh37", context_type = "default",
  #                      exome = False, minimum_signatures=1, maximum_signatures=10,
  #                      nmf_replicates=100, resample = True, batch_size=1, cpu=-1,
  #                      gpu=False, nmf_init="alexandrov-lab-custom", precision= "single",
  #                      matrix_normalization= "100X", seeds= "none",
  #                      min_nmf_iterations= 10000, max_nmf_iterations=1000000,
  #                      nmf_test_conv= 10000, nmf_tolerance= 1e-15, nnls_penalty=0.05,
  #                      get_all_signature_matrices= False)
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
        cpu = cores
      )
      sys$stdout$flush()
    },
    error = function(e) {
      print(e$message)
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
