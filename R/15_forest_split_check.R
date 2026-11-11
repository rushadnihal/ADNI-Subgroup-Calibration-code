# 15_forest_split_check.R ------------------------------------------------------
# Does record-level splitting inflate AUC for a flexible model? Script 14 found no inflation for
# Cox-B on stacked landmark rows (a low-capacity model cannot memorise subjects). A random survival
# forest can, so the same comparison is repeated for it:
#    RSF-default   500 trees, min node size 15  (the setting used in 07)
#    RSF-deep      500 trees, min node size 3   (high capacity, the case most likely to leak)
# Both use Tier C features (age, sex, education, APOE4, CDR-SB, ADAS-Cog13, RAVLT, FAQ, hippocampus,
# entorhinal, middle temporal, ventricles / ICV) plus a landmark index, fitted on the stacked landmark
# rows (2,732 rows, 1,194 subjects) with
#    (i)  the subject-level folds (all of a subject's rows in one fold), and
#    (ii) random record-level folds (the same subject's rows can fall on both sides).
# AUC (IPCW, KM censoring) is compared pooled at 2 and 3 years and by landmark at 2 years, with a
# bootstrap over subjects on the saved out-of-fold predictions (models are not refit).
# Cox-B results from script 14 are added to the figure for comparison when 14_split_auc_inflation.csv exists.
# Death is censoring. Counts and metrics only in the CSVs.
# Writes: output/audit/15_*.csv, 15_forest_split_report.txt, output/figures/15_forest_split.png (+pdf)
# Run:    source("R/15_forest_split_check.R")   (after 07 and 14; about 5 to 15 minutes)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
for (p in c("survival", "ranger"))
  if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is missing; run install.packages('", p, "').")
K <- 5; HORIZONS <- CFG$horizons_years
source(if (file.exists("R/utils_eval.R")) "R/utils_eval.R" else "utils_eval.R")

report <- file.path(CFG$dir_out, "audit", "15_forest_split_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Forest split check | run:", format(Sys.time()), "\n")

# ---- settings -------------------------------------------------------------
B   <- 500
FB  <- c("age", "male", "edu", "apoe4", "cdrsb", "adas13", "ravlt_immediate", "faq")
FC  <- c(FB, "hippo", "entorhinal", "midtemp", "ventricles")
FEATS <- c(FC, "lmi")
CONFIGS <- list(list(name = "RSF-default", mns = 15), list(name = "RSF-deep", mns = 3))
PAL <- c("#0072B2", "#D55E00", "#009E73")
# ---------------------------------------------------------------------------

tests <- data.frame(test = character(), result = character(), detail = character(), stringsAsFactors = FALSE)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s %s\n", if (ok) "PASS" else "FAIL", name, detail))
  tests <<- rbind(tests, data.frame(test = name, result = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE))
  invisible(ok)
}
rd <- function(x, k = 3) { nn <- vapply(x, is.numeric, logical(1)); x[nn] <- lapply(x[nn], round, k); x }

learner_rsf <- function(mns) function(dtr, dte, feats) {
  tr <- dtr[, c("time", "ev", feats)]; te <- dte[, feats, drop = FALSE]
  fit <- ranger::ranger(dependent.variable.name = "time", status.variable.name = "ev", data = tr,
                        num.trees = 500, min.node.size = mns, seed = CFG$seed, num.threads = 1)
  S <- predict(fit, data = te)$survival; ut <- fit$unique.death.times
  sapply(HORIZONS, function(h) { i <- sum(ut <= h); if (i == 0) rep(0, nrow(S)) else 1 - S[, i] })
}
run_cv_gen <- function(d, feats, learner) {
  n <- nrow(d); P <- matrix(NA_real_, n, length(HORIZONS))
  for (k in seq_len(K)) {
    tr <- d$fold != k; te <- !tr
    if (sum(d$ev[tr]) < 5 || sum(te) < 2) next
    num <- feats[sapply(feats, function(f) is.numeric(d[[f]]))]
    med <- sapply(num, function(f) median(d[[f]][tr], na.rm = TRUE))
    imp <- function(x) { for (f in num) x[[f]][is.na(x[[f]])] <- med[[f]]; x }
    res <- tryCatch(learner(imp(d[tr, ]), imp(d[te, ]), feats),
                    error = function(e) { cat("  note: fold", k, "failed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(res)) P[te, ] <- res
  }
  P
}
share_shared <- function(d) {
  out <- logical(nrow(d))
  for (k in 1:K) { te <- d$fold == k; out[te] <- d$RID[te] %in% d$RID[!te] }
  mean(out)
}

oc    <- readRDS(file.path(CFG$dir_out, "interim", "outcomes.rds"))
vis   <- readRDS(file.path(CFG$dir_out, "interim", "visits.rds"))
folds <- readRDS(file.path(CFG$dir_out, "interim", "folds.rds"))

d_all <- NULL; PP <- list(); share_s <- NA_real_; share_r <- NA_real_
leak_df <- NULL; diff_df <- NULL

section("0. Stacked landmark rows", {
  d0 <- make_d0(oc, vis, folds, "in_set")
  d_all <- d0[!is.na(d0$fold), ]
  d_all$lmi <- match(d_all$landmark, CFG$landmarks)
  cat("Rows:", nrow(d_all), "| subjects:", length(unique(d_all$RID)), "| events by 2 years:", sum(d_all$ev == 1 & d_all$time <= 2), "\n")
})

section("1. Out-of-fold forests with subject-level and record-level folds", {
  set.seed(CFG$seed + 1900)
  d_r <- d_all; d_r$fold <- sample(rep(1:K, length.out = nrow(d_all)))
  share_s <- share_shared(d_all); share_r <- share_shared(d_r)
  cat(sprintf("Share of test rows whose subject also has a training row: subject-level %.3f, record-level %.3f\n", share_s, share_r))
  for (cf in CONFIGS) {
    t0 <- Sys.time()
    lr <- learner_rsf(cf$mns)
    PP[[cf$name]] <- list(subject = run_cv_gen(d_all, FEATS, lr), record = run_cv_gen(d_r, FEATS, lr))
    cat(sprintf("  %s fitted (%.1f min)\n", cf$name, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
})

section("2. AUC by split type (bootstrap over subjects)", {
  set.seed(CFG$seed + 2000)
  n <- nrow(d_all)
  spec <- rbind(data.frame(scope = "pooled", h = c(2, 3), stringsAsFactors = FALSE),
                data.frame(scope = CFG$landmarks, h = 2, stringsAsFactors = FALSE))
  sel <- function(lmk, scope) if (scope == "pooled") rep(TRUE, length(lmk)) else lmk == scope
  auc_one <- function(tm, ev, P, lmk, scope, h) {
    ix <- which(sel(lmk, scope) & !is.na(P[, 1]))
    if (length(ix) < 20 || sum(ev[ix] == 1 & tm[ix] <= h) < 2) return(NA_real_)
    unname(met3(tm[ix], ev[ix], P[ix, match(h, HORIZONS)], h)["AUC"])
  }
  rows_by <- split(seq_len(n), d_all$RID); rids <- names(rows_by)
  lk <- list(); dfl <- list()
  for (cf in CONFIGS) {
    PM <- PP[[cf$name]]
    pt <- matrix(NA_real_, nrow(spec), 2)
    for (r in seq_len(nrow(spec))) for (v in 1:2)
      pt[r, v] <- auc_one(d_all$time, d_all$ev, PM[[v]], d_all$landmark, spec$scope[r], spec$h[r])
    bs <- array(NA_real_, c(B, nrow(spec), 2))
    for (b in 1:B) {
      ii <- unlist(rows_by[sample(rids, replace = TRUE)], use.names = FALSE)
      tm <- d_all$time[ii]; ev <- d_all$ev[ii]; lm <- d_all$landmark[ii]
      for (r in seq_len(nrow(spec))) for (v in 1:2)
        bs[b, r, v] <- auc_one(tm, ev, PM[[v]][ii, , drop = FALSE], lm, spec$scope[r], spec$h[r])
    }
    for (r in seq_len(nrow(spec))) {
      ix <- which(sel(d_all$landmark, spec$scope[r]))
      for (v in 1:2) {
        ci <- q95(bs[, r, v])
        lk[[length(lk) + 1]] <- data.frame(model = cf$name, split = c("subject", "record")[v], scope = spec$scope[r], horizon = spec$h[r],
          n = length(ix), events_by_h = sum(d_all$ev[ix] == 1 & d_all$time[ix] <= spec$h[r]), AUC = pt[r, v], lcl = ci[1], ucl = ci[2],
          stringsAsFactors = FALSE)
      }
      ci <- q95(bs[, r, 2] - bs[, r, 1])
      dfl[[length(dfl) + 1]] <- data.frame(model = cf$name, scope = spec$scope[r], horizon = spec$h[r],
        AUC_subject = pt[r, 1], AUC_record = pt[r, 2], inflation = pt[r, 2] - pt[r, 1], lcl = ci[1], ucl = ci[2], stringsAsFactors = FALSE)
    }
  }
  leak_df <- do.call(rbind, lk); diff_df <- do.call(rbind, dfl)
  save_csv(rd(leak_df), "15_forest_split_auc.csv"); save_csv(rd(diff_df), "15_forest_split_inflation.csv")
  cat("\nAUC inflation from record-level splitting (record minus subject):\n")
  show(rd(diff_df), 30)
  cat("\nCompare with Cox-B in 14_split_auc_inflation.csv (at most 0.002, intervals including 0).\n")
})

section("3. Figure", {
  fdir <- file.path(CFG$dir_out, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  cox <- tryCatch(as.data.frame(readr::read_csv(file.path(CFG$dir_out, "audit", "14_split_auc_inflation.csv"), show_col_types = FALSE)), error = function(e) NULL)
  if (!is.null(cox)) { cox$model <- "Cox-B"; cox <- cox[, c("model", "scope", "horizon", "AUC_subject", "AUC_record", "inflation", "lcl", "ucl")] }
  allm <- rbind(if (!is.null(cox)) cox, diff_df)
  mods <- unique(allm$model); scopes <- unique(paste(allm$scope, allm$horizon))
  fig <- function() {
    par(mar = c(7, 4.8, 3, 1))
    yr <- range(c(0, allm$lcl, allm$ucl), na.rm = TRUE)
    plot(NA, xlim = c(0.5, length(scopes) + 0.5), ylim = yr * c(1, 1.1), xaxt = "n", xlab = "",
         ylab = "AUC, record-level minus subject-level split", main = "Does record-level splitting inflate AUC?")
    title(xlab = "Scope (AUC horizon, years)", line = 5.5)
    axis(1, at = seq_along(scopes), labels = sub(" ", " (", paste0(scopes, ")")), las = 2, cex.axis = 0.75)
    abline(h = 0, lty = 2, col = "grey40")
    for (m in seq_along(mods)) for (s in seq_along(scopes)) {
      r <- allm[allm$model == mods[m] & paste(allm$scope, allm$horizon) == scopes[s], ]
      if (nrow(r) != 1 || !is.finite(r$lcl)) next
      x <- s + (m - (length(mods) + 1) / 2) * 0.22
      segments(x, r$lcl, x, r$ucl, col = PAL[m], lwd = 1.5); points(x, r$inflation, pch = 16, col = PAL[m])
    }
    legend("topleft", legend = mods, col = PAL[seq_along(mods)], pch = 16, bty = "n", cex = 0.8)
  }
  ok <- tryCatch({
    grDevices::png(file.path(fdir, "15_forest_split.png"), width = 8, height = 5, units = "in", res = 200); fig(); grDevices::dev.off()
    grDevices::pdf(file.path(fdir, "15_forest_split.pdf"), width = 8, height = 5); fig(); grDevices::dev.off(); TRUE
  }, error = function(e) { while (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE); cat("  figure not written:", conditionMessage(e), "\n"); FALSE })
  assign("fig_ok15", ok, envir = .GlobalEnv)
  cat("Figure:", if (ok) "15_forest_split.png and .pdf" else "none", "\n")
})

section("4. Checks", {
  ok <- !is.null(d_all) && length(PP) == length(CONFIGS) && !is.null(diff_df)
  check("data and forests were built", ok)
  if (ok) {
    check("subject-level folds never split a subject; record-level folds do", isTRUE(share_s == 0) && isTRUE(share_r > 0.3),
          sprintf("(%.3f and %.3f)", share_s, share_r))
    complete <- vapply(PP, function(x) !anyNA(x$subject) && !anyNA(x$record), logical(1))
    check("every row has an out-of-fold forest prediction", all(complete))
    a <- leak_df[leak_df$split == "subject" & leak_df$scope == "pooled" & leak_df$horizon == 2, ]
    check("subject-level forest AUC at 2 years is in a plausible range (0.75 to 0.95)", all(a$AUC > 0.75 & a$AUC < 0.95), sprintf("(%s)", paste(round(a$AUC, 3), collapse = ", ")))
    check("bootstrap intervals are finite for the pooled scopes", all(is.finite(diff_df$lcl[diff_df$scope == "pooled"]) & is.finite(diff_df$ucl[diff_df$scope == "pooled"])))
    ins <- diff_df$inflation >= diff_df$lcl - 0.02 & diff_df$inflation <= diff_df$ucl + 0.02
    check("inflation estimates sit inside their intervals (within 0.02)", all(ins, na.rm = TRUE))
    check("figure written", isTRUE(get0("fig_ok15", ifnotfound = FALSE)) && file.exists(file.path(CFG$dir_out, "figures", "15_forest_split.png")))
  }
  cat(sprintf("\n%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  save_csv(tests, "15_tests.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
