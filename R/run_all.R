# run_all.R ---------------------------------------------------------------
# Runs the whole pipeline, 01 to 15, each script in its own fresh R process (so one script's
# objects can never leak into the next), and writes a one-table summary of what happened.
#
# For every script it records: exit status, minutes, number of checks passed and failed (read
# from the script's own NN_tests.csv, only if the file was written during this run), and the
# number of sections that reported "SECTION FAILED" in the script's report.
#
# Usage (from the project folder, analysis1111, or its R/ subfolder):
#     source("R/run_all.R")                       # all scripts, stops at the first failure
#     run_all(from = 7, to = 15)                  # a range
#     run_all(stop_on_fail = FALSE)               # keep going after a failure
# or from a terminal:  Rscript R/run_all.R
#
# Total time: about 25 minutes (23 measured; 15, 07 and 14 are the longest; 08, 13 and 14 use 1000 bootstrap
# replicates). Nothing is run twice and nothing is skipped: the ADNIMERGE2 data must already be in
# "data and methods". Per-script console output goes to output/logs/run_NN.log; the summary goes to
# output/audit/run_all_summary.csv; package versions to output/logs/session_info.txt.
# -------------------------------------------------------------------------

run_all <- function(from = 1, to = 15, stop_on_fail = TRUE,
                    rscript = file.path(R.home("bin"), "Rscript")) {
  start_dir <- getwd()
  on.exit(setwd(start_dir), add = TRUE)
  root <- NULL
  for (p in c(start_dir, normalizePath(file.path(start_dir, ".."), mustWork = FALSE))) {
    if (dir.exists(file.path(p, "data and methods"))) { root <- normalizePath(p); break }
  }
  if (is.null(root)) stop("Cannot find the 'data and methods' folder. setwd() to the analysis1111 folder (or its R/ subfolder) first.")
  setwd(root)

  scripts <- sort(list.files("R", pattern = "^[0-9]{2}_.*\\.R$"))
  scripts <- scripts[!grepl("^00_", scripts)]
  num <- as.integer(substr(scripts, 1, 2))
  keep <- num >= from & num <= to
  scripts <- scripts[keep]; num <- num[keep]
  if (!length(scripts)) stop("No scripts found between ", from, " and ", to, " in ", file.path(root, "R"))

  audit <- file.path(root, "output", "audit"); logs <- file.path(root, "output", "logs")
  dir.create(audit, recursive = TRUE, showWarnings = FALSE); dir.create(logs, recursive = TRUE, showWarnings = FALSE)

  pk <- c("dplyr", "readr", "stringr", "survival", "glmnet", "ranger")
  ver <- vapply(pk, function(x) if (requireNamespace(x, quietly = TRUE)) as.character(utils::packageVersion(x)) else "not installed", character(1))
  writeLines(c(R.version.string, format(Sys.time()), "", paste(names(ver), ver)), file.path(logs, "session_info.txt"))
  if (any(ver == "not installed")) cat("Note: not installed:", paste(names(ver)[ver == "not installed"], collapse = ", "), "\n")

  rows <- list()
  cat("Running", length(scripts), "scripts from", root, "\n")
  for (i in seq_along(scripts)) {
    s <- scripts[i]; nn <- sprintf("%02d", num[i])
    t0 <- Sys.time()
    cat(sprintf("[%s] %s ... ", format(t0, "%H:%M:%S"), s)); flush.console()
    log <- file.path(logs, paste0("run_", nn, ".log"))
    status <- suppressWarnings(system2(rscript, c("--no-save", "--no-restore", shQuote(file.path("R", s))), stdout = log, stderr = log))
    status <- if (is.null(status)) NA_integer_ else as.integer(status)
    mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

    fresh <- function(f) file.exists(f) && as.numeric(file.mtime(f)) >= as.numeric(t0) - 1
    tf <- file.path(audit, paste0(nn, "_tests.csv"))
    n_pass <- NA_integer_; n_fail <- NA_integer_
    if (fresh(tf)) {
      tt <- tryCatch(utils::read.csv(tf, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(tt) && "result" %in% names(tt)) { n_pass <- sum(tt$result == "PASS"); n_fail <- sum(tt$result == "FAIL") }
    }
    sec_fail <- NA_integer_
    rp <- list.files(audit, pattern = sprintf("^%s_.*report\\.txt$", nn), full.names = TRUE)
    rp <- rp[vapply(rp, fresh, logical(1))]
    if (length(rp)) sec_fail <- sum(grepl("SECTION FAILED", readLines(rp[1], warn = FALSE), fixed = TRUE))
    ok <- isTRUE(status == 0L) && !isTRUE(n_fail > 0) && !isTRUE(sec_fail > 0)
    cat(sprintf("%s (%.1f min%s%s)\n", if (ok) "ok" else "PROBLEM", mins,
                if (!is.na(n_pass)) sprintf(", checks %d pass / %d fail", n_pass, n_fail) else "",
                if (isTRUE(sec_fail > 0)) sprintf(", %d section(s) failed", sec_fail) else ""))
    rows[[i]] <- data.frame(script = s, exit_status = status, minutes = round(mins, 1), checks_pass = n_pass, checks_fail = n_fail,
                            sections_failed = sec_fail, ok = ok, log = file.path("output", "logs", basename(log)), stringsAsFactors = FALSE)
    if (!ok && stop_on_fail) { cat("Stopping at", s, "- see", log, "and the script's report in output/audit.\n"); break }
  }
  res <- do.call(rbind, rows)
  utils::write.csv(res, file.path(audit, "run_all_summary.csv"), row.names = FALSE)
  cat(sprintf("\n%d of %d scripts ran cleanly; %.0f minutes in total. Summary: output/audit/run_all_summary.csv\n",
              sum(res$ok), nrow(res), sum(res$minutes)))
  invisible(res)
}

if (!exists("RUN_ALL_NO_AUTO")) run_all()
