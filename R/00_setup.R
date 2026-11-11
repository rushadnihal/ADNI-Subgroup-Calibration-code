# 00_setup.R ---------------------------------------------------------------
# Configuration and helpers. Every other script sources this automatically,
# so you never need to run it on its own.
#
# HOW TO RUN: set the working directory to the project folder (analysis1111),
# or to its R/ subfolder, then:  source("R/01_audit.R")
# -------------------------------------------------------------------------

find_root <- function() {
  cands <- c(getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  for (p in cands) {
    if (dir.exists(file.path(p, "data and methods"))) return(normalizePath(p))
  }
  stop("Cannot find the 'data and methods' folder. ",
       "setwd() to the analysis1111 folder (or its R/ subfolder) first.")
}

ROOT <- find_root()

CFG <- list(
  root           = ROOT,
  dir_data       = file.path(ROOT, "data and methods"),
  dir_rda        = file.path(ROOT, "data and methods", "ADNIMERGE2", "data"),
  dir_out        = file.path(ROOT, "output"),
  seed           = 20261111L,
  missing_codes  = c(-4, -5),                 # ADNI codes for "insufficient/missing"
  landmarks      = c("bl", "m12", "m24", "m36"),
  horizons_years = c(1, 2, 3)
)

for (d in c("audit", "interim", "logs")) {
  dir.create(file.path(CFG$dir_out, d), recursive = TRUE, showWarnings = FALSE)
}

# Packages: dplyr, readr, stringr only (more are added in later scripts).
req  <- c("dplyr", "readr", "stringr")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  stop("Install missing packages first:\n  install.packages(c(",
       paste0('"', miss, '"', collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})
options(dplyr.summarise.inform = FALSE)
set.seed(CFG$seed)

# ---- helpers --------------------------------------------------------------

# Safe date conversion (accepts Date, POSIXct, or ISO strings).
to_date <- function(x) suppressWarnings(as.Date(as.character(x)))

# Load one table from the ADNIMERGE2 package data (.rda); factors -> character.
load_rda <- function(name) {
  f <- file.path(CFG$dir_rda, paste0(name, ".rda"))
  if (!file.exists(f)) stop("Missing file: ", f)
  e   <- new.env()
  nm  <- load(f, envir = e)
  obj <- as.data.frame(e[[nm[1]]], stringsAsFactors = FALSE)
  obj[] <- lapply(obj, function(x) if (is.factor(x)) as.character(x) else x)
  obj
}

# Load one CSV from 'data and methods'.
read_adni_csv <- function(file) {
  f <- file.path(CFG$dir_data, file)
  if (!file.exists(f)) stop("Missing file: ", f)
  as.data.frame(readr::read_csv(f, guess_max = 1e6, show_col_types = FALSE,
                                progress = FALSE))
}

# Run one block; on error print the message and carry on (so one failure
# does not hide the rest of the report).
section <- function(title, expr) {
  cat("\n", strrep("=", 72), "\n", title, "\n", strrep("=", 72), "\n", sep = "")
  res <- tryCatch(expr, error = function(e) {
    cat("!! SECTION FAILED: ", conditionMessage(e), "\n", sep = "")
    NULL
  })
  invisible(res)
}

save_csv <- function(df, name, sub = "audit") {
  readr::write_csv(df, file.path(CFG$dir_out, sub, name), progress = FALSE)
  invisible(df)
}

show <- function(df, n = 60) {
  print(as.data.frame(utils::head(df, n)), row.names = FALSE)
  if (nrow(df) > n) cat("... ", nrow(df) - n, " more rows (see the CSV)\n", sep = "")
  invisible(df)
}

cat("Setup OK. Project root:", ROOT, "\n")
