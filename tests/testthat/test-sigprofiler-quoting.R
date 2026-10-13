# Tests for the hand-rolled `python -c "..."` fallback command of the
# SigProfiler interop bridge (sigprofiler_extract recall path).
#
# Reproduces the failure class of GitHub issue #470: on Windows, a path such
# as C:\Users\...\sigprofiler_input.txt was embedded un-escaped into the
# single-quoted python string literals of the recall command, so python
# interpreted `\Users` as a truncated \UXXXXXXXX unicode escape and the
# recall died with:
#   SyntaxError: (unicode error) 'unicodeescape' codec can't decode bytes ...
#     truncated \UXXXXXXXX escape
#
# Contract under test: quote_opt() must emit literals that decode back to the
# exact input value (backslash-escaped quotes handled explicitly), and every
# value that cannot be carried safely through the shell + python layers must
# be rejected fail-closed instead of producing a broken command.

# --- Verification helpers (test-local) -------------------------------------

# Minimal decoder for the body of a python single-quoted string literal.
# Raises a SyntaxError-style condition for exactly the escape hazards the
# recall command must never contain (issue #470 class) and returns the decoded
# text otherwise.
sp_py_decode <- function(body, what = "string literal") {
  chars <- strsplit(body, "")[[1L]]
  out <- character(0L)
  i <- 1L
  n <- length(chars)
  escape_chars <- c("U", "u", "x", "N")
  while (i <= n) {
    ch <- chars[i]
    if (identical(ch, "\\")) {
      if (i == n) {
        stop(sprintf("SyntaxError: %s ends with a dangling backslash", what),
             call. = FALSE)
      }
      nx <- chars[i + 1L]
      if (identical(nx, "\\")) {
        out <- c(out, "\\")
        i <- i + 2L
        next
      }
      if (length(match(nx, escape_chars)) > 0L) {
        stop(sprintf(paste0(
          "SyntaxError: (unicode error) 'unicodeescape' codec can't decode ",
          "bytes: truncated \\%sXXXXXXXX escape"
        ), nx), call. = FALSE)
      }
      stop(sprintf("SyntaxError: unsupported escape \\%s in generated %s",
                   nx, what), call. = FALSE)
    }
    out <- c(out, ch)
    i <- i + 1L
  }
  paste0(out, collapse = "")
}

# Strip the outer shell double quotes of the recall command and undo the
# backslash folding the given shell dialect performs inside double quotes.
#   windows -> cmd.exe passes content verbatim (percent expansion is rejected
#             outright by the builder, see quote_opt fail-closed tests).
#   posix   -> sh collapses \\ before $, `, ", and \ .
sp_shell_unwrap <- function(payload, dialect) {
  if (!grepl('^".*"$', payload)) {
    stop("payload is not wrapped in an outer double-quote pair")
  }
  body <- substr(payload, 2L, nchar(payload) - 1L)
  if (identical(dialect, "posix")) {
    chars <- strsplit(body, "")[[1L]]
    out <- character(0L)
    i <- 1L
    n <- length(chars)
    while (i <= n) {
      ch <- chars[i]
      if (identical(ch, "\\") && i < n &&
          length(match(chars[i + 1L], c("$", "`", "\"", "\\"))) > 0L) {
        out <- c(out, chars[i + 1L])
        i <- i + 2L
      } else {
        out <- c(out, ch)
        i <- i + 1L
      }
    }
    body <- paste0(out, collapse = "")
  }
  body
}

# Extract every single-quoted python literal from the (already shell-unwrapped)
# payload and decode each one.
sp_py_literals <- function(payload) {
  pieces <- strsplit(payload, "'")[[1L]]
  idx <- seq_along(pieces)
  literals <- pieces[(idx %% 2L) == 0L]
  vapply(literals, sp_py_decode, character(1L), USE.NAMES = FALSE)
}

# A real Windows temp path from the issue #470 traceback layout.
sp_win_tmp <- "C:\\Users\\jrc356\\AppData\\Local\\Temp\\RtmpuA1nnQ\\dir8e9c7fc35e75\\sigprofiler_input.txt"
sp_win_out <- "C:\\Users\\jrc356\\sigminer-work\\sigprofiler"

# --- quote_opt: escaping -----------------------------------------------------

test_that("quote_opt escapes Windows backslashes so the literal decodes to the input path (#470)", {
  out <- quote_opt(sp_win_tmp, os = "windows")
  expect_match(out, "^'.*'$")
  body <- substr(out, 2L, nchar(out) - 1L)
  # Raises the #470 SyntaxError class while the escaping is missing.
  decoded <- sp_py_decode(body, "input path literal")
  expect_equal(decoded, sp_win_tmp)
})

test_that("quote_opt keeps shell-benign values byte-identical", {
  expect_equal(
    quote_opt("GRCh37", opt = "reference_genome", os = "posix"),
    "reference_genome='GRCh37'"
  )
  expect_equal(
    quote_opt("matrix", os = "windows"),
    "'matrix'"
  )
})

test_that("quote_opt passes rm_quote scalars through un-quoted", {
  expect_equal(
    quote_opt(2, opt = "minimum_signatures", rm_quote = TRUE, os = "posix"),
    "minimum_signatures=2"
  )
  expect_equal(
    quote_opt("False", opt = "exome", rm_quote = TRUE, os = "posix"),
    "exome=False"
  )
})

# --- quote_opt: fail-closed rejections ---------------------------------------

test_that("quote_opt rejects values it cannot carry safely (fail-closed)", {
  # Embedded single quote would close the python literal early.
  expect_error(quote_opt("a'b", os = "windows"))
  # Embedded double quote would close the outer shell quote in both dialects.
  expect_error(quote_opt('a"b', os = "windows"))
  expect_error(quote_opt('a"b', os = "posix"))
  # Line breaks inject new commands into the emitted line.
  expect_error(quote_opt("a\nb", os = "windows"))
  expect_error(quote_opt("a\rb", os = "posix"))
  # cmd.exe expands %VAR% inside double quotes.
  expect_error(quote_opt("100%", os = "windows"))
  # sh expands $ and backtick inside double quotes.
  expect_error(quote_opt("a$b", os = "posix"))
  expect_error(quote_opt("a`b", os = "posix"))
  # A lone backslash cannot be carried through POSIX double quotes without
  # shell-dependent folding: reject instead of emitting an ambiguous command.
  expect_error(quote_opt("C:\\tmp\\x", os = "posix"))
})

# --- sigprofiler_fallback_command: the full command-assembly surface ---------

test_that("sigprofiler_fallback_command emits a decodable recall command on windows (#470 class)", {
  cmd <- sigprofiler_fallback_command(
    py_path = "C:/Users/jrc356/.conda/envs/sigprofiler_env/python.exe",
    output = sp_win_out,
    tmp_file = sp_win_tmp,
    genome_build = "GRCh37",
    sig_ranges = c(2, 5),
    nrun = 10L,
    is_exome = FALSE,
    init_method = "random",
    refit_plot = FALSE,
    cores = -1L,
    os = "windows"
  )
  expect_match(
    cmd,
    "^C:/Users/jrc356/\\.conda/envs/sigprofiler_env/python\\.exe -c \".*\"$"
  )
  expect_match(cmd, "sig\\.sigProfilerExtractor\\(")
  expect_match(cmd, "from SigProfilerExtractor import sigpro as sig;")

  payload <- sp_shell_unwrap(sub("^[^ ]+ -c ", "", cmd), "windows")
  decoded <- sp_py_literals(payload)
  # The two user paths must survive the round trip exactly.
  expect_equal(length(match(sp_win_out, decoded)), 1L)
  expect_equal(length(match(sp_win_tmp, decoded)), 1L)
  # rm_quote scalars appear un-quoted.
  expect_match(cmd, "minimum_signatures=2")
  expect_match(cmd, "maximum_signatures=5")
  expect_match(cmd, "nmf_replicates=10")
  expect_match(cmd, "exome=False")
  expect_match(cmd, "make_decomposition_plots=False")
  expect_match(cmd, "cpu=-1")
})

test_that("sigprofiler_fallback_command fails closed on shell-hostile py_path", {
  base <- list(
    output = "/tmp/sp-out", tmp_file = "/tmp/sp-in.txt",
    genome_build = "GRCh37", sig_ranges = c(2, 5), nrun = 10L,
    is_exome = FALSE, init_method = "random", refit_plot = FALSE,
    cores = -1L
  )
  # Space splits the command word.
  expect_error(do.call(sigprofiler_fallback_command,
                       c(list(py_path = "/usr/local/bin/my python", os = "unix"), base)))
  # Shell expansion in the command word.
  expect_error(do.call(sigprofiler_fallback_command,
                       c(list(py_path = "py$HOME", os = "unix"), base)))
  expect_error(do.call(sigprofiler_fallback_command,
                       c(list(py_path = "py%TEMP%exe", os = "windows"), base)))
  # Quotes anywhere are fatal in both dialects.
  expect_error(do.call(sigprofiler_fallback_command,
                       c(list(py_path = "py'exe", os = "windows"), base)))
  # Backslash in the command word: fine on cmd, folded by sh -> reject on posix.
  expect_error(do.call(sigprofiler_fallback_command,
                       c(list(py_path = "C:\\python.exe", os = "unix"), base)))
  # Empty command word.
  expect_error(do.call(sigprofiler_fallback_command,
                       c(list(py_path = "", os = "unix"), base)))
})

test_that("sigprofiler_fallback_command rejects backslash paths under posix", {
  expect_error(sigprofiler_fallback_command(
    py_path = "python",
    output = "C:\\sigminer-work\\sigprofiler",
    tmp_file = sp_win_tmp,
    genome_build = "GRCh37",
    sig_ranges = c(2, 5),
    nrun = 10L,
    is_exome = FALSE,
    init_method = "random",
    refit_plot = FALSE,
    cores = -1L,
    os = "unix"
  ))
})

test_that("no emitted payload carries a lone-backslash unicode escape hazard (#470 pin)", {
  cmd <- sigprofiler_fallback_command(
    py_path = "python",
    output = sp_win_out,
    tmp_file = sp_win_tmp,
    genome_build = "GRCh37",
    sig_ranges = c(2, 5),
    nrun = 10L,
    is_exome = TRUE,
    init_method = "nndsvd_min",
    refit_plot = TRUE,
    cores = 4L,
    os = "windows"
  )
  payload <- sp_shell_unwrap(sub("^[^ ]+ -c ", "", cmd), "windows")
  # Every emitted literal must decode cleanly.
  decoded <- sp_py_literals(payload)
  expect_equal(length(decoded), 6L) # matrix, output, tmp, ref genome, opp genome, init
  expect_equal(length(match(sp_win_out, decoded)), 1L)
})
