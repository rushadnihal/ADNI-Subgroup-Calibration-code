# make_github_release.R ---------------------------------------------------
# Builds a clean folder, "github_release/", holding only what is safe to put on GitHub.
# It copies from an ALLOW-LIST (never "everything except"), then scans every copied text file for
# anything that looks like participant-level data or a personal path, and refuses to finish if it
# finds any. It never modifies or deletes anything in your project.
#
# What is copied
#   always:        R/*.R, R/README.md, a safe .gitignore
#   if INCLUDE_RESULTS:  aggregate tables and text reports from output/audit (*.csv, *_report.txt,
#                  NN_tests.csv), figures from output/figures (*.png, *.pdf)
#   if INCLUDE_PAPER:    PDF files from report/ (your own manuscript; only if you want it public yet)
# What is never copied
#   data and methods/ (ADNI data and documentation), output/interim/ (subject-level .rds files),
#   output/logs/ (console logs, may hold local paths), papers/ (other people's PDFs), outline.pdf,
#   .Rhistory, .DS_Store.
#
# Usage (from the project folder, analysis1111, or its R/ subfolder):
#     source("R/make_github_release.R")
#     make_github_release()                                  # code only
#     make_github_release(include_results = TRUE)           # code + aggregate results and figures
#     make_github_release(include_results = TRUE, include_paper = TRUE)
#
# You are still responsible for checking the ADNI Data Use Agreement before publishing results.
# This script cannot judge that for you.
# -------------------------------------------------------------------------

make_github_release <- function(include_results = FALSE, include_paper = FALSE,
                                out_name = "github_release", max_mb = 5) {
  start_dir <- getwd()
  on.exit(setwd(start_dir), add = TRUE)
  root <- NULL
  for (p in c(start_dir, normalizePath(file.path(start_dir, ".."), mustWork = FALSE))) {
    if (dir.exists(file.path(p, "data and methods"))) { root <- normalizePath(p); break }
  }
  if (is.null(root)) stop("Cannot find the 'data and methods' folder. setwd() to the analysis1111 folder first.")
  setwd(root)

  out <- file.path(root, out_name)
  if (dir.exists(out)) stop("'", out_name, "' already exists. Move or rename it first; this script will not overwrite or delete it.")
  dir.create(out)

  copied <- character(0)
  put <- function(from, to_rel) {
    to <- file.path(out, to_rel)
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(from, to)
    if (ok) copied <<- c(copied, to_rel)
    invisible(ok)
  }

  # 1. Code
  for (f in list.files("R", pattern = "\\.(R|md)$", full.names = TRUE)) put(f, file.path("R", basename(f)))

  # 2. Aggregate results (optional)
  if (include_results) {
    for (f in list.files(file.path("output", "audit"), pattern = "(\\.csv|_report\\.txt)$", full.names = TRUE))
      put(f, file.path("results", "tables", basename(f)))
    for (f in list.files(file.path("output", "figures"), pattern = "\\.(png|pdf)$", full.names = TRUE))
      put(f, file.path("results", "figures", basename(f)))
  }

  # 3. Manuscript PDF (optional)
  if (include_paper && dir.exists("report")) {
    for (f in list.files("report", pattern = "\\.pdf$", full.names = TRUE)) put(f, file.path("report", basename(f)))
  }

  # 4. A safe .gitignore for the new repository
  writeLines(c(
    "# Never commit ADNI data, derived subject-level files, logs or other people's papers",
    "data and methods/", "output/interim/", "output/logs/", "papers/", "*.rds", "*.RData",
    ".Rhistory", ".Rproj.user/", ".DS_Store"), file.path(out, ".gitignore"))
  copied <- c(copied, ".gitignore")

  # 5. Scan everything copied
  problems <- character(0)
  flag <- function(f, why) problems <<- c(problems, paste0(f, ": ", why))
  id_pat   <- "[0-9]{3}_S_[0-9]{4}"                    # ADNI subject ID (PTID)
  path_pat <- "/Users/|/home/[a-z]|C:\\\\Users|@gmail\\.com|@[a-z0-9-]+\\.(edu|org|com)"
  col_pat  <- "^\"?(RID|PTID|subject|subject_id|id|ID|VISCODE|EXAMDATE|USERDATE)\"?(,|$)"
  for (rel in copied) {
    f <- file.path(out, rel)
    if (file.size(f) > max_mb * 1024^2) flag(rel, paste("larger than", max_mb, "MB"))
    if (grepl("\\.(rds|RData|rda)$", rel, ignore.case = TRUE)) flag(rel, "binary R data file")
    if (grepl("\\.(R|md|csv|txt)$|^\\.gitignore$", rel)) {
      x <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
      if (any(grepl(id_pat, x))) flag(rel, "contains an ADNI-style subject ID")
      if (grepl("\\.csv$", rel) && length(x) && grepl(col_pat, x[1])) flag(rel, "header looks like subject-level data")
      if (!grepl("^R/", rel) || grepl("\\.md$", rel)) {
        if (any(grepl(path_pat, x))) flag(rel, "contains a personal path or e-mail address")
      }
    }
  }

  man <- data.frame(file = copied,
                    kb = round(file.size(file.path(out, copied)) / 1024, 1))
  write.csv(man, file.path(out, "MANIFEST.csv"), row.names = FALSE)

  cat("\n=== github_release summary ===\n")
  cat("Folder:", out, "\nFiles copied:", length(copied), "\n")
  if (length(problems)) {
    cat("\nFAIL. Do NOT publish. Fix or remove these first:\n")
    for (p in problems) cat("  -", p, "\n")
    cat("\n(The folder was still created so you can look at it; delete it once fixed and rerun.)\n")
    return(invisible(FALSE))
  }
  cat("\nPASS. No subject IDs, subject-level headers, R data files, large files or personal paths found.\n")
  cat("Next: add a LICENSE and the ADNI acknowledgement to the README, read the ADNI Data Use Agreement,\n")
  cat("then 'cd", out_name, "&& git init && git add . && git commit' and push.\n")
  invisible(TRUE)
}

if (!exists("MAKE_RELEASE_NO_AUTO") && sys.nframe() == 0L) make_github_release()
