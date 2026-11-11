# 10_plasma_competing_risk.R -------------------------------------------------
# Two sensitivity analyses.
#
# A. Tier D (plasma biomarkers) on the subset of landmark visits that have a
#    plasma draw within the look-back window (365 days). Plasma exists only for a
#    minority of visits, so this part (1) reports who has plasma, by subgroup and
#    within protocol era, (2) compares Tier B with and without plasma on exactly
#    the same rows, with paired bootstrap AUC gains, and (3) runs a
#    shuffled-plasma negative control. Landmarks bl and m12 only: later landmarks
#    have too few plasma draws (see 02).
#
# B. Death as a competing event. Scripts 04 to 09 treat death as censoring. Here
#    death within DEATH_WINDOW_DAYS of the last visit (set in 03) is its own
#    event. Compares 1 - Kaplan-Meier with the Aalen-Johansen cumulative
#    incidence, and compares Cox-B (death censored), a cause-specific Cox, and a
#    Fine-Gray model on discrimination and calibration against the
#    Aalen-Johansen observed risk. AUC treats deaths before the horizon as
#    controls, weighted by inverse probability of censoring (not by death).
#
# Reporting rules: subgroup levels with fewer than MIN_N subjects are not
# analysed. Folds are the subject-level folds saved by 04.
#
# Writes: output/audit/10_*.csv, 10_plasma_competing_risk_report.txt
#         output/figures/10_*.png (aggregate results only)
# Run:    source("R/10_plasma_competing_risk.R")   (after 07; a few minutes)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
for (p in c("survival", "glmnet"))
  if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is missing; run install.packages('", p, "').")
K <- 5; HORIZONS <- CFG$horizons_years
source(if (file.exists("R/utils_eval.R")) "R/utils_eval.R" else "utils_eval.R")

report <- file.path(CFG$dir_out, "audit", "10_plasma_competing_risk_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Plasma tier and competing risks | run:", format(Sys.time()), "\n")

# ---- settings -------------------------------------------------------------
B      <- 500      # bootstrap replicates, part A
B_CR   <- 300      # bootstrap replicates, part B
MIN_N  <- 30
FT0    <- c("age", "apoe4", "cdrsb")
FB     <- c("age", "male", "edu", "apoe4", "cdrsb", "adas13", "ravlt_immediate", "faq")
PL_RAW <- c("ptau217", "ab42_ab40", "nfl", "gfap")
PL     <- c("lptau217", "ab42_ab40", "lnfl", "lgfap")     # log of the right-skewed markers
PL_LM  <- c("bl", "m12")
VARS   <- c("sex", "race_bw", "hisp", "edu_grp", "apoe4_grp", "age_grp", "origprot")
# ---------------------------------------------------------------------------

tests <- data.frame(test = character(), result = character(), detail = character(), stringsAsFactors = FALSE)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s %s\n", if (ok) "PASS" else "FAIL", name, detail))
  tests <<- rbind(tests, data.frame(test = name, result = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE))
  invisible(ok)
}
rd <- function(x) { nn <- vapply(x, is.numeric, logical(1)); x[nn] <- lapply(x[nn], round, 3); x }

# ---- learners and helpers (own copies, identical to 07 where shared) -------------
breslow <- function(time, ev, lp, hs) {
  o <- order(time); time <- time[o]; ev <- ev[o]; r <- exp(lp[o])
  rs <- rev(cumsum(rev(r)))
  cum <- cumsum(ifelse(ev == 1, 1 / rs, 0))
  sapply(hs, function(h) { i <- sum(time <= h); if (i == 0) 0 else cum[i] })
}
learn_glmnet <- function(dtr, dte, feats) {
  mm <- function(x) model.matrix(~ . - 1, data = x[, feats, drop = FALSE])
  xtr <- mm(dtr); xte <- mm(dte)
  ytr <- cbind(time = dtr$time, status = dtr$ev)
  set.seed(CFG$seed)
  cvf <- glmnet::cv.glmnet(xtr, ytr, family = "cox", alpha = 0.5, nfolds = 5)
  lam <- cvf$lambda.min
  lp_tr <- as.numeric(predict(cvf, newx = xtr, s = lam, type = "link"))
  lp_te <- as.numeric(predict(cvf, newx = xte, s = lam, type = "link"))
  m <- mean(lp_tr); lp_tr <- lp_tr - m; lp_te <- lp_te - m
  H0 <- breslow(dtr$time, dtr$ev, lp_tr, HORIZONS)
  sapply(H0, function(h0) 1 - exp(-h0 * exp(lp_te)))
}
# Fine-Gray: cumulative incidence of dementia (state 1) with death (state 2) competing.
fg_fit_predict <- function(dtr, dte, feats) {
  df <- data.frame(time = dtr$time_cr,
                   st = factor(dtr$status_cr, levels = 0:2, labels = c("cens", "dem", "death")),
                   dtr[, feats, drop = FALSE])
  fg  <- survival::finegray(survival::Surv(time, st) ~ ., data = df, etype = "dem")
  fit <- survival::coxph(stats::as.formula(paste("survival::Surv(fgstart, fgstop, fgstatus) ~", paste(feats, collapse = " + "))),
                         data = fg, weights = fgwt)
  sf  <- survival::survfit(fit, newdata = dte[, feats, drop = FALSE])
  sapply(HORIZONS, function(h) 1 - as.numeric(summary(sf, times = h, extend = TRUE)$surv))
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
# Aalen-Johansen cumulative incidence of dementia at hs (falls back to 1 - KM when nobody died).
aj_risk <- function(time, status, hs) {
  if (sum(status == 1) == 0) return(rep(0, length(hs)))
  if (!any(status == 2))
    return(1 - summary(survival::survfit(survival::Surv(time, status == 1) ~ 1), times = hs, extend = TRUE)$surv)
  sf <- survival::survfit(survival::Surv(time, factor(status, levels = 0:2, labels = c("cens", "dem", "death"))) ~ 1)
  as.numeric(summary(sf, times = hs, extend = TRUE)$pstate[, match("dem", sf$states)])
}
# Competing-risk AUC at h: cases = dementia by h; controls = still free at h, or dead before h and
# before dementia. Censoring weights come from a Kaplan-Meier of loss to follow-up (status 0 only).
auc_cr <- function(time, status, p, h) {
  ok <- !is.na(p); time <- time[ok]; status <- status[ok]; p <- p[ok]
  G <- cens_G(time, as.integer(status != 0))
  case <- status == 1 & time <= h; ctrl <- time > h; dth <- status == 2 & time <= h
  if (!any(case) || !any(ctrl | dth)) return(NA_real_)
  wc <- 1 / G(time[case]); nk <- ctrl | dth
  wk <- ifelse(ctrl[nk], 1 / G(h), 1 / G(time[nk]))
  pc <- p[case]; pk <- p[nk]
  sum(outer(wc, wk) * ((outer(pc, pk, ">")) + 0.5 * outer(pc, pk, "=="))) / (sum(wc) * sum(wk))
}
make_specs <- function(dl) {
  n <- nrow(dl); specs <- list(list(var = "all", level = "all", idx = seq_len(n)))
  for (v in VARS) {
    lv <- as.character(dl[[v]])
    for (u in setdiff(unique(lv), c(NA, "Unknown", "missing")))
      if (sum(lv == u, na.rm = TRUE) >= MIN_N) specs[[length(specs) + 1]] <- list(var = v, level = u, idx = which(lv == u))
  }
  specs
}
lg <- function(x) { m <- suppressWarnings(min(x[x > 0], na.rm = TRUE)); log(pmax(x, m / 2)) }

# ---- data -----------------------------------------------------------------------
oc    <- readRDS(file.path(CFG$dir_out, "interim", "outcomes.rds"))
vis   <- readRDS(file.path(CFG$dir_out, "interim", "visits.rds"))
folds <- readRDS(file.path(CFG$dir_out, "interim", "folds.rds"))
o07   <- readRDS(file.path(CFG$dir_out, "interim", "oof_predictions.rds"))
stopifnot(all(PL_RAW %in% names(vis)), identical(as.numeric(o07$horizons), as.numeric(HORIZONS)))

d0 <- make_d0(oc, vis, folds, "in_set")
pv <- vis[, c("RID", "EXAMDATE", PL_RAW)]
mm <- match(paste(d0$RID, d0$t0), paste(pv$RID, pv$EXAMDATE))
for (v in PL_RAW) d0[[v]] <- pv[[v]][mm]
d0$lptau217 <- lg(d0$ptau217); d0$lnfl <- lg(d0$nfl); d0$lgfap <- lg(d0$gfap)
d0$has_pl <- stats::complete.cases(d0[, PL_RAW])
m2 <- match(paste(d0$RID, d0$landmark), paste(oc$RID, oc$landmark))
d0$time_cr <- oc$time[m2]; d0$status_cr <- oc$status[m2]
stopifnot(!anyNA(d0$time_cr), !anyNA(d0$status_cr))
aligned <- identical(paste(d0$RID, d0$landmark), paste(o07$d$RID, o07$d$landmark))
cat("Analysis rows:", nrow(d0), "| aligned with saved predictions:", aligned, "\n")

# ============================ PART A: plasma ==========================================
section("A1. Who has a plasma draw within 365 days of the landmark visit", {
  cov_rows <- list(); era_rows <- list()
  for (l in PL_LM) {
    dl <- d0[d0$landmark == l, ]
    for (s in make_specs(dl)) {
      i <- s$idx
      cov_rows[[length(cov_rows) + 1]] <- data.frame(landmark = l, variable = s$var, level = s$level, n = length(i),
        n_plasma = sum(dl$has_pl[i]), share_plasma = mean(dl$has_pl[i]), stringsAsFactors = FALSE)
    }
    for (e in levels(dl$era)) for (rc in c("Black", "White")) {
      i <- which(dl$era == e & dl$race_bw == rc)
      era_rows[[length(era_rows) + 1]] <- data.frame(landmark = l, era = e, race = rc, n = length(i),
        n_plasma = sum(dl$has_pl[i]), share_plasma = if (length(i)) mean(dl$has_pl[i]) else NA_real_, stringsAsFactors = FALSE)
    }
  }
  cov_df <- rd(do.call(rbind, cov_rows)); era_df <- rd(do.call(rbind, era_rows))
  save_csv(cov_df, "10_plasma_coverage_by_subgroup.csv"); save_csv(era_df, "10_plasma_coverage_race_within_era.csv")
  assign("cov_df", cov_df, envir = .GlobalEnv)
  show(cov_df, 80)
  cat("\nBlack vs White share with plasma, within each protocol era:\n")
  show(era_df, 40)
  for (l in PL_LM) {
    dl <- d0[d0$landmark == l & d0$race_bw %in% c("Black", "White"), ]
    tb <- table(has_pl = dl$has_pl, race = dl$race_bw, era = dl$era)
    r <- tryCatch(stats::mantelhaen.test(tb), error = function(e) NULL)
    if (!is.null(r)) cat(sprintf("%s: Mantel-Haenszel odds ratio of having plasma, Black vs White, stratified by era: %.2f (%.2f to %.2f), p = %.3f\n",
                                 l, 1 / r$estimate, 1 / r$conf.int[2], 1 / r$conf.int[1], r$p.value))
  }
  cat("\nRows with plasma at the later landmarks (not analysed):\n")
  show(data.frame(landmark = setdiff(CFG$landmarks, PL_LM),
                  n = sapply(setdiff(CFG$landmarks, PL_LM), function(l) sum(d0$landmark == l)),
                  n_plasma = sapply(setdiff(CFG$landmarks, PL_LM), function(l) sum(d0$landmark == l & d0$has_pl))))
})

PA <- list(); DA <- list(); mainA <- list(); deltaA <- list()
section("A2. Tier B with and without plasma, on the same rows (out-of-fold)", {
  for (l in PL_LM) {
    dl <- d0[d0$landmark == l & d0$has_pl, ]; rownames(dl) <- NULL; n <- nrow(dl)
    P <- list()
    P[["Cox-T0"]]          <- run_cv(dl, FT0)$P
    P[["Cox-T0+pTau217"]]  <- run_cv(dl, c(FT0, "lptau217"))$P
    P[["Cox-B"]]           <- run_cv(dl, FB)$P
    P[["Cox-B+plasma"]]    <- run_cv(dl, c(FB, PL))$P
    P[["GLMNET-B+plasma"]] <- run_cv_gen(dl, c(FB, PL), learn_glmnet)
    ix <- match(dl$RID, o07$d$RID[o07$d$landmark == l])   # row within this landmark's saved matrix
    P[["Cox-B (full-cohort fit)"]] <- o07$P[[l]][["Cox-B"]][ix, , drop = FALSE]
    # negative control: shuffle the plasma rows jointly
    set.seed(CFG$seed + 1000 + match(l, CFG$landmarks))
    ds <- dl; pm <- sample.int(n)
    for (v in PL) ds[[v]] <- dl[[v]][pm]
    P[["Cox-B+shuffled plasma"]] <- run_cv(ds, c(FB, PL))$P
    PA[[l]] <- P; DA[[l]] <- dl
    mn <- names(P)
    set.seed(CFG$seed + 1100 + match(l, CFG$landmarks))
    bt <- array(NA_real_, c(B, length(mn), length(HORIZONS), 3),
                dimnames = list(NULL, mn, as.character(HORIZONS), c("AUC", "sBrier", "OE")))
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      for (m in mn) for (j in seq_along(HORIZONS))
        bt[b, m, j, ] <- met3(dl$time[idx], dl$ev[idx], P[[m]][idx, j], HORIZONS[j])
    }
    for (m in mn) for (j in seq_along(HORIZONS)) {
      h <- HORIZONS[j]; pt <- met3(dl$time, dl$ev, P[[m]][, j], h)
      mainA[[length(mainA) + 1]] <- data.frame(landmark = l, model = m, horizon = h, n = n,
        events_by_h = sum(dl$ev == 1 & dl$time <= h), AUC = pt[["AUC"]],
        AUC_lcl = q95(bt[, m, j, "AUC"])[1], AUC_ucl = q95(bt[, m, j, "AUC"])[2],
        sBrier = pt[["sBrier"]], OE = pt[["OE"]], OE_lcl = q95(bt[, m, j, "OE"])[1], OE_ucl = q95(bt[, m, j, "OE"])[2],
        stringsAsFactors = FALSE)
      if (m != "Cox-B") {
        d_pt <- pt[["AUC"]] - met3(dl$time, dl$ev, P[["Cox-B"]][, j], h)[["AUC"]]
        cd <- q95(bt[, m, j, "AUC"] - bt[, "Cox-B", j, "AUC"])
        deltaA[[length(deltaA) + 1]] <- data.frame(landmark = l, model = m, horizon = h, events_by_h = sum(dl$ev == 1 & dl$time <= h),
          delta_AUC_vs_CoxB = d_pt, lcl = cd[1], ucl = cd[2], stringsAsFactors = FALSE)
      }
    }
    cat("  ", l, ": n =", n, "| events by 1/2/3 years:", paste(sapply(HORIZONS, function(h) sum(dl$ev == 1 & dl$time <= h)), collapse = "/"), "\n")
  }
  mainA_df <- rd(do.call(rbind, mainA)); deltaA_df <- rd(do.call(rbind, deltaA))
  save_csv(mainA_df, "10_plasma_metrics.csv"); save_csv(deltaA_df, "10_plasma_auc_gain.csv")
  assign("mainA_df", mainA_df, envir = .GlobalEnv); assign("deltaA_df", deltaA_df, envir = .GlobalEnv)
  cat("\nAUC, scaled Brier and O/E on the plasma subset, 2-year horizon:\n")
  show(mainA_df[mainA_df$horizon == 2, c("landmark", "model", "n", "events_by_h", "AUC", "AUC_lcl", "AUC_ucl", "sBrier", "OE", "OE_lcl", "OE_ucl")], 40)
  cat("\nAUC change relative to Cox-B trained on the subset (paired bootstrap):\n")
  show(deltaA_df[deltaA_df$horizon %in% c(1, 2, 3), ], 60)
})

# ============================ PART B: competing risks =================================
section("B1. Death as a competing event: counts and observed risk", {
  rows <- list(); cnt <- list()
  for (l in CFG$landmarks) {
    dl <- d0[d0$landmark == l, ]
    for (h in HORIZONS) cnt[[length(cnt) + 1]] <- data.frame(landmark = l, horizon = h, n = nrow(dl),
      dementia_by_h = sum(dl$status_cr == 1 & dl$time_cr <= h), death_first_by_h = sum(dl$status_cr == 2 & dl$time_cr <= h))
    for (s in make_specs(dl)) {
      i <- s$idx
      aj <- aj_risk(dl$time_cr[i], dl$status_cr[i], HORIZONS)
      for (j in seq_along(HORIZONS)) {
        h <- HORIZONS[j]; km <- km_risk(dl$time[i], dl$ev[i], h)
        rows[[length(rows) + 1]] <- data.frame(landmark = l, variable = s$var, level = s$level, horizon = h, n = length(i),
          dementia_by_h = sum(dl$status_cr[i] == 1 & dl$time_cr[i] <= h), death_first_by_h = sum(dl$status_cr[i] == 2 & dl$time_cr[i] <= h),
          risk_1_minus_KM = km, risk_AJ = aj[j], difference_pp = 100 * (km - aj[j]), stringsAsFactors = FALSE)
      }
    }
  }
  cnt_df <- do.call(rbind, cnt); cr_df <- rd(do.call(rbind, rows))
  save_csv(cr_df, "10_cr_risk_comparison.csv"); assign("cr_df", cr_df, envir = .GlobalEnv)
  cat("Dementia and death-first counts by horizon (death counts only inside the", 365, "day window set in 03):\n")
  show(cnt_df, 20)
  cat("\nObserved risk of dementia: 1 - Kaplan-Meier (death censored) vs Aalen-Johansen (death competing), whole cohort:\n")
  show(cr_df[cr_df$variable == "all", c("landmark", "horizon", "dementia_by_h", "death_first_by_h", "risk_1_minus_KM", "risk_AJ", "difference_pp")], 20)
  big <- cr_df[cr_df$variable != "all" & cr_df$horizon == 2 & abs(cr_df$difference_pp) >= 1, ]
  cat("\nSubgroups where the 2-year risks differ by 1 percentage point or more (point estimates):\n")
  if (nrow(big)) show(big[, c("landmark", "variable", "level", "n", "dementia_by_h", "death_first_by_h", "risk_1_minus_KM", "risk_AJ", "difference_pp")], 40) else cat("  none\n")
})

CRM <- list(); mainB <- list(); deltaB <- list()
section("B2. Cox-B with death censored, cause-specific Cox, and Fine-Gray", {
  for (l in CFG$landmarks) {
    dl <- d0[d0$landmark == l, ]; rownames(dl) <- NULL; n <- nrow(dl)
    ix <- match(dl$RID, o07$d$RID[o07$d$landmark == l])   # row within this landmark's saved matrix
    dcs <- dl; dcs$time <- dl$time_cr; dcs$ev <- as.integer(dl$status_cr == 1)
    P <- list("Cox-B (death censored)" = o07$P[[l]][["Cox-B"]][ix, , drop = FALSE],
              "Cox-B cause-specific"   = run_cv(dcs, FB)$P,
              "Fine-Gray B"            = run_cv_gen(dl, FB, fg_fit_predict))
    CRM[[l]] <- P; mn <- names(P)
    set.seed(CFG$seed + 1200 + match(l, CFG$landmarks))
    bt <- array(NA_real_, c(B_CR, length(mn), length(HORIZONS), 2), dimnames = list(NULL, mn, as.character(HORIZONS), c("AUC", "OE")))
    for (b in seq_len(B_CR)) {
      idx <- sample.int(n, n, replace = TRUE); aj <- aj_risk(dl$time_cr[idx], dl$status_cr[idx], HORIZONS)
      for (m in mn) for (j in seq_along(HORIZONS)) {
        pj <- P[[m]][idx, j]
        bt[b, m, j, ] <- c(auc_cr(dl$time_cr[idx], dl$status_cr[idx], pj, HORIZONS[j]), aj[j] / mean(pj, na.rm = TRUE))
      }
    }
    aj0 <- aj_risk(dl$time_cr, dl$status_cr, HORIZONS)
    for (m in mn) for (j in seq_along(HORIZONS)) {
      h <- HORIZONS[j]; pj <- P[[m]][, j]
      auc <- auc_cr(dl$time_cr, dl$status_cr, pj, h); oe <- aj0[j] / mean(pj, na.rm = TRUE)
      mainB[[length(mainB) + 1]] <- data.frame(landmark = l, model = m, horizon = h, n = n,
        dementia_by_h = sum(dl$status_cr == 1 & dl$time_cr <= h), AJ_observed = aj0[j], mean_predicted = mean(pj, na.rm = TRUE),
        AUC_cr = auc, AUC_lcl = q95(bt[, m, j, "AUC"])[1], AUC_ucl = q95(bt[, m, j, "AUC"])[2],
        OE_vs_AJ = oe, OE_lcl = q95(bt[, m, j, "OE"])[1], OE_ucl = q95(bt[, m, j, "OE"])[2], stringsAsFactors = FALSE)
      if (m != mn[1]) {
        base <- auc_cr(dl$time_cr, dl$status_cr, P[[mn[1]]][, j], h)
        cd <- q95(bt[, m, j, "AUC"] - bt[, mn[1], j, "AUC"])
        deltaB[[length(deltaB) + 1]] <- data.frame(landmark = l, model = m, horizon = h, delta_AUC_vs_censored = auc - base, lcl = cd[1], ucl = cd[2], stringsAsFactors = FALSE)
      }
    }
    cat("  ", l, ": done\n")
  }
  mainB_df <- rd(do.call(rbind, mainB)); deltaB_df <- rd(do.call(rbind, deltaB))
  save_csv(mainB_df, "10_cr_models.csv"); save_csv(deltaB_df, "10_cr_auc_gain.csv")
  assign("mainB_df", mainB_df, envir = .GlobalEnv); assign("deltaB_df", deltaB_df, envir = .GlobalEnv)
  cat("\nObserved risk is Aalen-Johansen. 2-year horizon:\n")
  show(mainB_df[mainB_df$horizon == 2, c("landmark", "model", "dementia_by_h", "AJ_observed", "mean_predicted", "AUC_cr", "AUC_lcl", "AUC_ucl", "OE_vs_AJ", "OE_lcl", "OE_ucl")], 40)
  cat("\nAUC difference from Cox-B with death censored (paired bootstrap):\n")
  show(deltaB_df, 60)
})

# ============================ figures ================================================
section("Figures", {
  fdir <- file.path(CFG$dir_out, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  ok1 <- tryCatch({
    f <- file.path(fdir, "10_plasma_auc.png")
    grDevices::png(f, width = 1600, height = 700, res = 150)
    par(mfrow = c(1, 2), mar = c(4, 11, 3, 1))
    for (l in PL_LM) {
      d <- mainA_df[mainA_df$landmark == l & mainA_df$horizon == 2, ]; d <- d[nrow(d):1, ]; yy <- seq_len(nrow(d))
      plot(NA, xlim = c(0.4, 1), ylim = c(0.5, nrow(d) + 0.5), yaxt = "n", xlab = "2-year AUC (out-of-fold)", ylab = "",
           main = paste0(l, " (n = ", d$n[1], ", ", d$events_by_h[1], " events)"), cex.main = 0.9)
      axis(2, at = yy, labels = d$model, las = 1, cex.axis = 0.65); abline(v = 0.5, lty = 3, col = "grey60")
      segments(d$AUC_lcl, yy, d$AUC_ucl, yy, col = "#3b6ea5", lwd = 1.5); points(d$AUC, yy, pch = 19, col = "#3b6ea5", cex = 0.8)
    }
    grDevices::dev.off(); TRUE
  }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); cat("  figure 1 not written:", conditionMessage(e), "\n"); FALSE })
  ok2 <- tryCatch({
    f <- file.path(fdir, "10_km_vs_aj_risk.png")
    grDevices::png(f, width = 1600, height = 700, res = 150)
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    d <- cr_df[cr_df$variable == "all", ]
    plot(NA, xlim = c(0.5, 4.5), ylim = range(c(d$risk_1_minus_KM, d$risk_AJ)) * c(0.9, 1.1), xaxt = "n", xlab = "", ylab = "Observed 3-year risk", main = "Risk: 1 - KM (solid) vs Aalen-Johansen (dashed)", cex.main = 0.85)
    axis(1, at = 1:4, labels = CFG$landmarks)
    for (h in HORIZONS) { x <- d[d$horizon == h, ]; x <- x[order(match(x$landmark, CFG$landmarks)), ]
      lines(1:4, x$risk_1_minus_KM, pch = 19, type = "o", col = "#3b6ea5"); lines(1:4, x$risk_AJ, pch = 1, type = "o", lty = 2, col = "#c0504d") }
    mtext("Lines: horizons 1, 2, 3 years (low to high)", side = 1, line = 2.5, cex = 0.6)
    plot(NA, xlim = c(0.5, 4.5), ylim = range(c(0, d$difference_pp)) * 1.1 + c(0, 0.1), xaxt = "n", xlab = "", ylab = "1 - KM minus AJ (percentage points)", main = "Overstatement of risk when death is censored", cex.main = 0.85)
    axis(1, at = 1:4, labels = CFG$landmarks); abline(h = 0, col = "grey50")
    for (h in HORIZONS) { x <- d[d$horizon == h, ]; x <- x[order(match(x$landmark, CFG$landmarks)), ]; lines(1:4, x$difference_pp, type = "o", pch = 19, col = grDevices::hcl.colors(3, "Dark 3")[match(h, HORIZONS)]) }
    legend("topleft", legend = paste(HORIZONS, "year"), col = grDevices::hcl.colors(3, "Dark 3"), pch = 19, lty = 1, bty = "n", cex = 0.7)
    grDevices::dev.off(); TRUE
  }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); cat("  figure 2 not written:", conditionMessage(e), "\n"); FALSE })
  cat("Figures written to output/figures/ :", c("10_plasma_auc.png", "10_km_vs_aj_risk.png")[c(ok1, ok2)], "\n")
})

# ============================ checks ================================================
section("Checks", {
  check("analysis rows line up with the saved predictions", aligned)
  r8 <- tryCatch(readr::read_csv(file.path(CFG$dir_out, "audit", "08_oe_over_time.csv"), show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(r8)) {
    a <- r8[r8$variable == "all" & r8$model == "Cox-B", ]
    ok <- sapply(seq_len(nrow(a)), function(k) { d1 <- d0[d0$landmark == a$landmark[k], ]; abs(km_risk(d1$time, d1$ev, a$horizon[k]) - a$km_observed[k]) < 0.002 })
    check("whole-cohort 1 - KM risk reproduces 08 for every landmark and horizon (within 0.002)", all(ok), sprintf("(%d rows)", nrow(a)))
  } else check("08 results found for the regression check", FALSE)
  check("every plasma-subset row is found in the saved predictions for its landmark",
        all(sapply(PL_LM, function(l) !anyNA(match(DA[[l]]$RID, o07$d$RID[o07$d$landmark == l])))))
  check("plasma subset uses only rows with all four markers present",
        all(sapply(PL_LM, function(l) !anyNA(DA[[l]][, PL_RAW]) && nrow(DA[[l]]) == sum(d0$landmark == l & d0$has_pl))))
  # competing-risk machinery on synthetic data with known cumulative incidence
  set.seed(5); n <- 4000; x <- rnorm(n); l1 <- 0.2 * exp(0.7 * x); l2 <- 0.1
  t1 <- rexp(n, l1); t2 <- rexp(n, l2); tc <- runif(n, 0.5, 6); tt <- pmin(t1, t2, tc)
  st <- ifelse(tc < pmin(t1, t2), 0L, ifelse(t1 < t2, 1L, 2L))
  true2 <- mean(l1 / (l1 + l2) * (1 - exp(-(l1 + l2) * 2))); j2 <- match(2, HORIZONS)
  aj2 <- aj_risk(tt, st, 2); km2 <- 1 - summary(survival::survfit(survival::Surv(tt, st == 1) ~ 1), times = 2, extend = TRUE)$surv
  check("synthetic: Aalen-Johansen recovers the true cumulative incidence and 1 - KM overstates it",
        abs(aj2 - true2) < 0.03 && km2 > aj2, sprintf("(true %.3f, AJ %.3f, 1-KM %.3f)", true2, aj2, km2))
  sd <- data.frame(time_cr = tt, status_cr = st, x = x)
  fgp <- tryCatch(fg_fit_predict(sd, sd, "x")[, j2], error = function(e) NA_real_)
  check("synthetic: Fine-Gray mean predicted risk matches the true cumulative incidence (within 0.03)",
        isTRUE(abs(mean(fgp) - true2) < 0.03), sprintf("(true %.3f, predicted %.3f)", true2, mean(fgp)))
  y0 <- as.integer(st == 1 & tt <= 2)
  check("synthetic: competing-risk AUC is 1 for a perfect predictor, near 0.5 for noise",
        { a1 <- auc_cr(tt, st, as.numeric(y0), 2); a2 <- auc_cr(tt, st, runif(n), 2); isTRUE(a1 > 0.99 && abs(a2 - 0.5) < 0.06) },
        sprintf("(perfect %.3f, noise %.3f)", auc_cr(tt, st, as.numeric(y0), 2), auc_cr(tt, st, runif(n), 2)))
  chk <- all(sapply(CFG$landmarks, function(l) !anyNA(CRM[[l]][["Fine-Gray B"]]) ))
  check("Fine-Gray produced predictions for every subject at every landmark", chk)
  sh <- deltaA_df[deltaA_df$model == "Cox-B+shuffled plasma" & deltaA_df$landmark == "bl" & deltaA_df$horizon == 2, ]
  check("negative control: shuffled plasma does not raise baseline 2-year AUC by more than 0.03", nrow(sh) == 1 && sh$delta_AUC_vs_CoxB < 0.03,
        sprintf("(change %.3f)", if (nrow(sh)) sh$delta_AUC_vs_CoxB else NA))
  cat(sprintf("\n%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  save_csv(tests, "10_tests.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
