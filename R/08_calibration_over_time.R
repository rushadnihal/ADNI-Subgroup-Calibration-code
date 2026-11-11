# 08_calibration_over_time.R -----------------------------------------------
# Subgroup calibration tracked over landmark time, from the out-of-fold
# predictions saved by 07 (no models are refit).
#
# For every landmark (bl, m12, m24, m36), horizon (1, 2, 3 years), model and
# subgroup level: O/E = Kaplan-Meier observed risk / mean predicted risk, with
# subject-level bootstrap intervals. Also bootstrap ratios of O/E between pairs
# of subgroup levels (Black vs White, Hispanic vs not, and so on): a ratio whose
# interval excludes 1 means the model is calibrated differently for the two groups.
#
# Reporting rules: levels with fewer than MIN_N subjects are listed but not
# analysed; levels with fewer than MIN_EV events by the horizon get a point
# estimate but no interval and are flagged. Death is treated as censoring.
#
# Writes: output/audit/08_*.csv, 08_calibration_over_time_report.txt
#         output/figures/08_*.png (counts and ratios only)
# Run:    source("R/08_calibration_over_time.R")   (after 07; about 5 to 10 minutes)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
if (!requireNamespace("survival", quietly = TRUE))
  stop("Package 'survival' is missing. It ships with R; run install.packages('survival').")
K <- 5; HORIZONS <- CFG$horizons_years
source(if (file.exists("R/utils_eval.R")) "R/utils_eval.R" else "utils_eval.R")

report <- file.path(CFG$dir_out, "audit", "08_calibration_over_time_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Calibration over landmark time | run:", format(Sys.time()), "\n")

# ---- settings -------------------------------------------------------------
B      <- 1000
MIN_N  <- 30
MIN_EV <- 5
VARS   <- c("sex", "race_bw", "hisp", "edu_grp", "apoe4_grp", "age_grp", "origprot")
CONTRASTS <- list(c("race_bw", "Black", "White"), c("hisp", "Hispanic", "Not Hispanic"),
                  c("sex", "Female", "Male"), c("edu_grp", "<=12y", ">12y"),
                  c("apoe4_grp", "1+", "0"), c("age_grp", "75+", "<65"),
                  c("origprot", "ADNI3", "ADNI1"), c("origprot", "ADNI4", "ADNI1"))
FIG_MODEL <- "Cox-B"; FIG_H <- 2
# ---------------------------------------------------------------------------

tests <- data.frame(test = character(), result = character(), detail = character(), stringsAsFactors = FALSE)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s %s\n", if (ok) "PASS" else "FAIL", name, detail))
  tests <<- rbind(tests, data.frame(test = name, result = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE))
  invisible(ok)
}
km_risk_h <- function(time, ev) {
  if (sum(ev) == 0) return(rep(0, length(HORIZONS)))
  1 - summary(survival::survfit(survival::Surv(time, ev) ~ 1), times = HORIZONS, extend = TRUE)$surv
}

dat <- readRDS(file.path(CFG$dir_out, "interim", "oof_predictions.rds"))
D <- dat$d; Pm <- dat$P
stopifnot(identical(as.numeric(dat$horizons), as.numeric(HORIZONS)))

res <- list(); con <- list(); skipped <- list()
section("1. Subgroup O/E at every landmark, horizon and model", {
  for (l in CFG$landmarks) {
    dl <- D[D$landmark == l, ]; n <- nrow(dl); Pl <- Pm[[l]]; mods <- names(Pl)
    # subgroup specifications (levels with enough subjects)
    specs <- list(list(var = "all", level = "all", lv = rep("all", n)))
    for (v in VARS) {
      lv <- as.character(dl[[v]])
      for (u in setdiff(unique(lv), c(NA, "Unknown", "missing"))) {
        k <- sum(lv == u, na.rm = TRUE)
        if (k >= MIN_N) specs[[length(specs) + 1]] <- list(var = v, level = u, lv = lv)
        else skipped[[length(skipped) + 1]] <- data.frame(landmark = l, variable = v, level = u, n = k)
      }
    }
    Q <- length(specs)
    idx_of <- lapply(specs, function(s) which(s$lv == s$level))
    # point estimates
    obs <- t(sapply(idx_of, function(s) km_risk_h(dl$time[s], dl$ev[s])))              # Q x H
    pred <- lapply(mods, function(m) t(sapply(idx_of, function(s) colMeans(Pl[[m]][s, , drop = FALSE], na.rm = TRUE))))
    names(pred) <- mods
    evn <- t(sapply(idx_of, function(s) sapply(HORIZONS, function(h) sum(dl$ev[s] == 1 & dl$time[s] <= h))))
    # bootstrap
    oe_b <- array(NA_real_, c(B, Q, length(mods), length(HORIZONS)))
    set.seed(CFG$seed + 300 + match(l, CFG$landmarks))
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      for (q in seq_len(Q)) {
        s <- idx[which(specs[[q]]$lv[idx] == specs[[q]]$level)]
        if (length(s) < 5) next
        ob <- km_risk_h(dl$time[s], dl$ev[s])
        for (mi in seq_along(mods)) oe_b[b, q, mi, ] <- ob / colMeans(Pl[[mods[mi]]][s, , drop = FALSE], na.rm = TRUE)
      }
      if (b %% 250 == 0) cat("  ", l, ": bootstrap", b, "of", B, "\n")
    }
    for (q in seq_len(Q)) for (mi in seq_along(mods)) for (j in seq_along(HORIZONS)) {
      oe <- obs[q, j] / pred[[mods[mi]]][q, j]
      low <- evn[q, j] < MIN_EV
      ci <- if (low) c(NA_real_, NA_real_) else q95(oe_b[, q, mi, j])
      res[[length(res) + 1]] <- data.frame(landmark = l, model = mods[mi], horizon = HORIZONS[j],
        variable = specs[[q]]$var, level = specs[[q]]$level, n = length(idx_of[[q]]), events_by_h = evn[q, j],
        km_observed = obs[q, j], mean_pred = pred[[mods[mi]]][q, j], O_over_E = oe, lcl = ci[1], ucl = ci[2],
        reported = !low, stringsAsFactors = FALSE)
    }
    # contrasts between subgroup levels
    for (ct in CONTRASTS) {
      q1 <- which(sapply(specs, function(s) s$var == ct[1] && s$level == ct[2]))
      q2 <- which(sapply(specs, function(s) s$var == ct[1] && s$level == ct[3]))
      if (length(q1) != 1 || length(q2) != 1) next
      for (mi in seq_along(mods)) for (j in seq_along(HORIZONS)) {
        low <- evn[q1, j] < MIN_EV || evn[q2, j] < MIN_EV
        pt <- (obs[q1, j] / pred[[mods[mi]]][q1, j]) / (obs[q2, j] / pred[[mods[mi]]][q2, j])
        ci <- if (low) c(NA_real_, NA_real_) else q95(oe_b[, q1, mi, j] / oe_b[, q2, mi, j])
        con[[length(con) + 1]] <- data.frame(landmark = l, model = mods[mi], horizon = HORIZONS[j],
          variable = ct[1], group = ct[2], reference = ct[3], n_group = length(idx_of[[q1]]), n_reference = length(idx_of[[q2]]),
          events_group = evn[q1, j], events_reference = evn[q2, j], ratio_of_OE = pt, lcl = ci[1], ucl = ci[2],
          excludes_1 = if (low) NA else (ci[1] > 1 || ci[2] < 1), stringsAsFactors = FALSE)
      }
    }
  }
  rd <- function(x) { nn <- vapply(x, is.numeric, logical(1)); x[nn] <- lapply(x[nn], round, 3); x }
  res_df <- rd(do.call(rbind, res)); con_df <- rd(do.call(rbind, con))
  sk_df <- if (length(skipped)) do.call(rbind, skipped) else data.frame(landmark = character(), variable = character(), level = character(), n = integer())
  save_csv(res_df, "08_oe_over_time.csv"); save_csv(con_df, "08_oe_contrasts.csv"); save_csv(sk_df, "08_levels_too_small.csv")
  assign("res_df", res_df, envir = .GlobalEnv); assign("con_df", con_df, envir = .GlobalEnv)
  cat("\nWhole cohort, O/E by landmark and horizon (all models in 08_oe_over_time.csv):\n")
  w <- res_df[res_df$variable == "all" & res_df$model %in% c("Cox-B", "GLMNET-C", "RSF-C"),
              c("landmark", "model", "horizon", "events_by_h", "O_over_E", "lcl", "ucl")]
  show(w, 50)
  cat("\nLevels with fewer than", MIN_N, "subjects (not analysed):", nrow(sk_df), "- see 08_levels_too_small.csv\n")
})

section("2. Ratio of O/E between subgroup levels (model Cox-B, 2-year horizon)", {
  x <- con_df[con_df$model == "Cox-B" & con_df$horizon == 2,
              c("landmark", "variable", "group", "reference", "events_group", "events_reference", "ratio_of_OE", "lcl", "ucl", "excludes_1")]
  show(x, 80)
  r <- con_df[!is.na(con_df$excludes_1), ]
  cat(sprintf("\nAcross all models, horizons and landmarks: %d of %d contrasts with enough events have an interval that excludes 1.\n",
              sum(r$excludes_1), nrow(r)))
  if (any(r$excludes_1)) {
    cat("They are:\n")
    show(r[r$excludes_1, c("landmark", "model", "horizon", "variable", "group", "reference", "ratio_of_OE", "lcl", "ucl")], 80)
  }
  cat("\nWith about 40 contrasts per model and horizon, expect a few intervals to exclude 1 by chance (95% intervals, no multiplicity correction).\n")
})

section("3. Figures", {
  fdir <- file.path(CFG$dir_out, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  clip <- function(x) pmin(pmax(x, 0.25), 4)
  panel <- function(df, v, ttl) {
    d <- df[df$variable == v & df$reported, ]
    lvs <- unique(d$level); cols <- grDevices::hcl.colors(max(3, length(lvs)), "Dark 3")[seq_along(lvs)]
    plot(NA, xlim = c(0.5, 4.5), ylim = c(0.25, 4), log = "y", xaxt = "n", xlab = "", ylab = "O/E (log scale)", main = ttl, cex.main = 0.9)
    axis(1, at = 1:4, labels = CFG$landmarks); abline(h = 1, lty = 2, col = "grey40")
    off <- if (length(lvs) > 1) seq(-0.25, 0.25, length.out = length(lvs)) else 0
    for (i in seq_along(lvs)) {
      x <- d[d$level == lvs[i], ]; x <- x[order(match(x$landmark, CFG$landmarks)), ]
      xx <- match(x$landmark, CFG$landmarks) + off[i]
      segments(xx, clip(x$lcl), xx, clip(x$ucl), col = cols[i], lwd = 1.5)
      lines(xx, clip(x$O_over_E), col = cols[i]); points(xx, clip(x$O_over_E), pch = 19, col = cols[i], cex = 0.8)
    }
    if (length(lvs)) legend("topright", legend = lvs, col = cols, pch = 19, bty = "n", cex = 0.65)
  }
  ok1 <- tryCatch({
    f <- file.path(fdir, "08_oe_over_time_by_subgroup.png")
    grDevices::png(f, width = 2000, height = 1100, res = 150)
    par(mfrow = c(2, 4), mar = c(3, 4, 3, 1))
    d <- res_df[res_df$model == FIG_MODEL & res_df$horizon == FIG_H, ]
    for (v in c("all", VARS)) panel(d, v, paste0(v, " (", FIG_MODEL, ", ", FIG_H, "-year)"))
    grDevices::dev.off(); TRUE
  }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); cat("  figure 1 not written:", conditionMessage(e), "\n"); FALSE })
  ok2 <- tryCatch({
    f <- file.path(fdir, "08_overall_oe_by_model.png")
    grDevices::png(f, width = 2000, height = 700, res = 150)
    par(mfrow = c(1, length(HORIZONS)), mar = c(3, 4, 3, 1))
    for (h in HORIZONS) {
      d <- res_df[res_df$variable == "all" & res_df$horizon == h, ]
      mods <- unique(d$model); cols <- grDevices::hcl.colors(length(mods), "Dark 3")
      plot(NA, xlim = c(0.5, 4.5), ylim = c(0.6, 1.5), xaxt = "n", xlab = "", ylab = "O/E", main = paste0("Whole cohort, ", h, "-year horizon"), cex.main = 0.9)
      axis(1, at = 1:4, labels = CFG$landmarks); abline(h = 1, lty = 2, col = "grey40")
      for (i in seq_along(mods)) {
        x <- d[d$model == mods[i], ]; x <- x[order(match(x$landmark, CFG$landmarks)), ]
        xx <- match(x$landmark, CFG$landmarks) + (i - (length(mods) + 1) / 2) * 0.08
        segments(xx, pmax(x$lcl, 0.6), xx, pmin(x$ucl, 1.5), col = cols[i]); lines(xx, x$O_over_E, col = cols[i]); points(xx, x$O_over_E, pch = 19, col = cols[i], cex = 0.7)
      }
      if (h == HORIZONS[1]) legend("topleft", legend = mods, col = cols, pch = 19, bty = "n", cex = 0.55)
    }
    grDevices::dev.off(); TRUE
  }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); cat("  figure 2 not written:", conditionMessage(e), "\n"); FALSE })
  cat("Figures written to output/figures/ :", c("08_oe_over_time_by_subgroup.png", "08_overall_oe_by_model.png")[c(ok1, ok2)], "\n")
})

section("4. Checks", {
  ref7 <- tryCatch(readr::read_csv(file.path(CFG$dir_out, "audit", "07_model_metrics.csv"), show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(ref7)) {
    a <- res_df[res_df$variable == "all", c("landmark", "model", "horizon", "O_over_E", "events_by_h")]
    m <- merge(a, as.data.frame(ref7[, c("landmark", "model", "horizon", "O_over_E", "events_by_h")]),
               by = c("landmark", "model", "horizon"), suffixes = c("_08", "_07"))
    check("whole-cohort O/E reproduces 07 for every landmark, model and horizon (within 0.005)",
          nrow(m) == nrow(ref7) && all(abs(m$O_over_E_08 - m$O_over_E_07) < 0.005), sprintf("(%d rows, largest difference %.4f)", nrow(m), max(abs(m$O_over_E_08 - m$O_over_E_07))))
    check("event counts match 07 for every landmark, model and horizon", all(m$events_by_h_08 == m$events_by_h_07))
  } else check("07 results found for the regression check", FALSE)
  check("saved predictions line up with the outcome rows at every landmark",
        all(sapply(CFG$landmarks, function(l) nrow(D[D$landmark == l, ]) == nrow(Pm[[l]][[1]]))))
  rr <- res_df[res_df$reported & !is.na(res_df$lcl), ]
  check("bootstrap intervals contain their point estimates (within 5% of O/E)",
        mean(rr$O_over_E >= rr$lcl * 0.95 & rr$O_over_E <= rr$ucl * 1.05) > 0.99, sprintf("(%.1f%% of %d rows)", 100 * mean(rr$O_over_E >= rr$lcl * 0.95 & rr$O_over_E <= rr$ucl * 1.05), nrow(rr)))
  check("contrast ratios equal the ratio of the two O/E values",
        { x <- con_df[con_df$model == "Cox-B" & con_df$horizon == 2 & con_df$landmark == "bl" & con_df$variable == "race_bw", ]
          a1 <- res_df[res_df$model == "Cox-B" & res_df$horizon == 2 & res_df$landmark == "bl" & res_df$variable == "race_bw", ]
          nrow(x) == 1 && abs(x$ratio_of_OE - a1$O_over_E[a1$level == "Black"] / a1$O_over_E[a1$level == "White"]) < 0.01 })
  cat(sprintf("\n%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  save_csv(tests, "08_tests.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
