# Shared path and workflow helpers for the HiCBZIP release scripts.

hicbzip_this_file <- function() {
  frames <- sys.frames()
  files <- vapply(frames, function(x) {
    if (!is.null(x$ofile)) normalizePath(x$ofile, winslash = "/", mustWork = FALSE) else NA_character_
  }, character(1))
  files <- files[!is.na(files)]
  if (length(files) == 0) return(NA_character_)
  files[[length(files)]]
}

HICBZIP_COMMON_DIR <- dirname(hicbzip_this_file())
HICBZIP_ROOT <- normalizePath(file.path(HICBZIP_COMMON_DIR, ".."), winslash = "/", mustWork = FALSE)

path_here <- function(...) {
  file.path(HICBZIP_ROOT, ...)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

require_files <- function(paths, label = "required file") {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "Missing ", label, ":\n- ",
      paste(normalizePath(missing, winslash = "/", mustWork = FALSE), collapse = "\n- "),
      call. = FALSE
    )
  }
  invisible(paths)
}

require_dirs <- function(paths, label = "required directory") {
  missing <- paths[!dir.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "Missing ", label, ":\n- ",
      paste(normalizePath(missing, winslash = "/", mustWork = FALSE), collapse = "\n- "),
      call. = FALSE
    )
  }
  invisible(paths)
}

require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(pkgs)
}

source_hicbzip_core <- function(include_score_helpers = FALSE, include_embedding_helpers = FALSE) {
  require_files(c(
    path_here("HiCBZIP", "analysis_helpers.R"),
    path_here("HiCBZIP", "BHZIP_EB.R")
  ), label = "HiCBZIP core source")

  source(path_here("HiCBZIP", "analysis_helpers.R"))
  source(path_here("HiCBZIP", "BHZIP_EB.R"))

  if (include_score_helpers) {
    require_files(path_here("_common", "score_io_helpers.R"), label = "SCORE helper source")
    source(path_here("_common", "score_io_helpers.R"))
  }
  if (include_embedding_helpers) {
    require_files(path_here("_common", "embedding_helpers.R"), label = "embedding helper source")
    source(path_here("_common", "embedding_helpers.R"))
  }

  invisible(TRUE)
}

render_rmd <- function(input, output_dir = NULL, ...) {
  require_packages("rmarkdown")
  require_files(input, label = "R Markdown workflow")
  if (!is.null(output_dir)) ensure_dir(output_dir)
  rmarkdown::render(input = input, output_dir = output_dir, envir = new.env(parent = globalenv()), ...)
}
