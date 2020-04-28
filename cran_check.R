rhub::check_for_cran(path = ".", email = NULL, check_args = "--as-cran",
                     env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false",
                                  `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "true",
                                  R_COMPILE_AND_INSTALL_PACKAGES = "true"))
