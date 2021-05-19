#' Get Curated Reference Signature Database
#'
#' Reference mutational signatures and their aetiologies,
#' mainly obtained from COSMIC database
#' (SigProfiler results) and cleaned before saving into
#' **sigminer** package. You can obtain:
#' - COSMIC legacy SBS signatures.
#' - COSMIC v3 SBS signatures.
#' - COSMIC v3 DBS signatures.
#' - COSMIC v3 ID (indel) signatures.
#' - SBS and RS (rearrangement) signatures from Nik lab 2020 Nature Cancer paper.
#' - RS signatures from BRCA560 and USARC cohorts.
#' @inheritParams get_sig_similarity
#'
#' @return a `list`.
#' @export
#' @seealso [get_sig_similarity], [sig_fit] and [show_cosmic_sig_profile].
#'
#' @examples
#' s1 <- get_sig_db()
#' s2 <- get_sig_db("SBS")
#' s3 <- get_sig_db("DBS")
#' s4 <- get_sig_db("DBS_mm10")
#' s5 <- get_sig_db("SBS_Nik_lab")
#' s6 <- get_sig_db("ID")
#' s7 <- get_sig_db("RS_BRCA560")
#' s8 <- get_sig_db("RS_USARC")
#' s9 <- get_sig_db("RS_Nik_lab")
#' s1
#' s2
#' s3
#' s4
#' s5
#' s6
#' s7
#' s8
#' s9
#' @testexamples
#' expect_is(s1, "list")
#' expect_is(s2, "list")
#' expect_is(s3, "list")
#' expect_is(s4, "list")
#' expect_is(s5, "list")
#' expect_is(s6, "list")
#' expect_is(s7, "list")
#' expect_is(s8, "list")
#' expect_is(s9, "list")
get_sig_db <- function(sig_db = "legacy") {
  db_file <- switch(sig_db,
    legacy = system.file("extdata", "legacy_signatures.RDs",
      package = "maftools", mustWork = TRUE
    ),
    SBS = system.file("extdata", "SBS_signatures.rds",
      package = "sigminer", mustWork = TRUE
    ),
    SBS_Nik_lab = system.file("extdata", "SBS_signatures_Nik_lab.rds",
      package = "sigminer", mustWork = TRUE
    ),
    SBS_Nik_lab_Organ = system.file("extdata", "SBS_signatures_Nik_lab_Organ.rds",
      package = "sigminer", mustWork = TRUE
    ),
    RS_Nik_lab = system.file("extdata", "RS_signatures_Nik_lab.rds",
      package = "sigminer", mustWork = TRUE
    ),
    RS_Nik_lab_Organ = system.file("extdata", "RS_signatures_Nik_lab_Organ.rds",
      package = "sigminer", mustWork = TRUE
    ),
    RS_BRCA560 = system.file("extdata", "RS_signatures_BRCA560.rds",
      package = "sigminer", mustWork = TRUE
    ),
    RS_USARC = system.file("extdata", "RS_signatures_USARC.rds",
      package = "sigminer", mustWork = TRUE
    ),
    CNS_USARC = system.file("extdata", "CNS_signatures_USARC.rds",
      package = "sigminer", mustWork = TRUE
    ),
    DBS = system.file("extdata", "DBS_signatures.rds",
      package = "sigminer", mustWork = TRUE
    ),
    ID = system.file("extdata", "ID_signatures.rds",
      package = "sigminer", mustWork = TRUE
    ),
    TSB = system.file("extdata", "TSB_signatures.rds",
      package = "sigminer", mustWork = TRUE
    ),
    {
      if (startsWith(sig_db, "SBS")) {
        system.file("extdata", "SBS_signatures_list.rds",
          package = "sigminer", mustWork = TRUE
        )
      } else if (startsWith(sig_db, "DBS")) {
        system.file("extdata", "DBS_signatures_list.rds",
          package = "sigminer", mustWork = TRUE
        )
      } else if (startsWith(sig_db, "latest_")) {
        latest_version <- "v3.2"
        # v3.2_SBS_GRCh37.txt
        sig_db2 <- sub("latest_", "", sig_db)

        file_dir <- switch(sig_db,
          latest_SBS_GRCh37 = "452",
          latest_DBS_GRCh37 = "446",
          latest_ID_GRCh37 = "451",
          latest_SBS_GRCh38 = "453",
          latest_DBS_GRCh38 = "447",
          latest_SBS_mm9 = "454",
          latest_DBS_mm9 = "448",
          latest_SBS_mm10 = "455",
          latest_DBS_mm10 = "449",
          latest_SBS_rn6 = "456",
          latest_DBS_rn6 = "450",
          stop("Bad input!")
        )

        data_url <- file.path(
          "https://cancer.sanger.ac.uk/signatures/documents", file_dir,
          paste0("COSMIC_", latest_version, "_", sig_db2, ".txt")
        )
        db_file <- file.path(
          system.file("extdata", package = "sigminer"),
          paste0("COSMIC_", latest_version, "_", sig_db2, ".rds")
        )
        # 读入后需要处理下
        if (!file.exists(db_file)) {
          message("The data is not available in local, obtain it from COSMIC: https://cancer.sanger.ac.uk/signatures/downloads/")
          data <- tryCatch(
            data.table::fread(data_url, data.table = FALSE) %>%
              tibble::column_to_rownames("Type") %>%
              as.matrix(),
            error = function(e) {
              message(e$message)
              stop("The data cannot be downloaded, please check your internet or contact the developer.")
            }
          )
          message("Transforming and saving to ", db_file)
          adata <- list(
            db = data,
            aetiology = data.frame(
              V1 = colnames(data),
              V2 = "See COSMIC https://cancer.sanger.ac.uk/signatures/"
            ) %>% tibble::column_to_rownames("V1") %>% setNames("aetiology"),
            date = latest_version
          )
          saveRDS(adata, file = db_file)
        }
        db_file
      } else {
        stop("Bad input for option 'sig_db'!")
      }
    }
  )
  sigs_db <- readRDS(file = db_file)
  if (!"db" %in% names(sigs_db)) {
    sigs_db$db <- sigs_db[[sub(".*_", "", sig_db)]]
    sigs_db <- sigs_db[c("db", "aetiology", "date")]
  }
  ## Make sure column-sum is 1, i.e. normalized
  sigs_db$db <- apply(sigs_db$db, 2, function(x) x / sum(x))
  sigs_db
}
