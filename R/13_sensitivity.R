# 13_sensitivity.R -------------------------------------------------------------
# Two sensitivity analyses that the results draft still flags.
#
# A. CONFIRMED-DEMENTIA OUTCOME. 22 of the 201 baseline conversions within 2 years revert at
#    the next visit. Here the event is the first dementia diagnosis that is confirmed (the next
#    visit is also dementia, or it is the subject's last visit), the event time is that visit
#    date, and death is still treated as censoring. Three scenarios, same folds, Cox on Tier B:
#      S1  main outcome,      model trained on the main outcome        (the reference)
#      S2  confirmed outcome, model trained on the main outcome        (does the headline survive a stricter outcome?)
#      S3  confirmed outcome, model trained on the confirmed outcome   (and does retraining change it?)
#    For each: AUC, scaled Brier and O/E by landmark and horizon, subgroup O/E by protocol era,
#    race and ethnicity, and the ADNI3-vs-ADNI1 and Black-vs-White ratios of O/E, with subject-level
#    bootstrap intervals (one resample shared by all scenarios, so differences are paired).
#
# B. CENSORING BOUNDS FOR COX-B. 06 bounded how far loss to follow-up could explain the era and race
#    findings using the age/APOE4/CDR-SB model. Repeated here with the Cox-B predictions saved by 07:
#      B1  observed risk if none / all of the lost subjects had converted, and the tipping fraction
#      B2  covariate-adjusted IPCW observed risk versus Kaplan-Meier
#    Groups: protocol era (ADNI1, ADNI-GO/2, ADNI3/4), each protocol, race, ethnicity.
#
# Death is censoring throughout. Counts and metrics only in the CSVs.
# Writes: output/audit/13_*.csv, 13_sensitivity_report.txt, output/figures/13_era_oe_by_outcome.png (+pdf)
# Run:    source("R/13_sensitivity.R")   (after 07; about 15 to 25 minutes)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
if (!requireNamespace("survival", quietly = TRUE))
  stop("Package 'survival' is missing. It ships with R; run install.packages('survival').")
K <- 5; HORIZONS <- CFG$horizons_years
source(if (file.exists("R/utils_eval.R")) "R/utils_eval.R" else "utils_eval.R")

report <- file.path(CFG$dir_out, "audit", "13_sensitivity_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Sensitivity analyses | run:", format(Sys.time()), "\n")

# ---- settings -------------------------------------------------------------
B        <- 1000
MIN_N    <- 30
MIN_EV   <- 5
DAYS     <- 365.25
FB       <- c("age", "male", "edu", "apoe4", "cdrsb", "adas13", "ravlt_immediate", "faq")
LM_ALL   <- CFG$landmarks
LM_SUB   <- c("bl", "m12")
HS_SUB   <- c(2, 3)
SUBV     <- c("origprot", "race_bw", "hisp")
LM_B     <- c("bl", "m12")
HS_B     <- c(2, 3)
GRP_B    <- c("all", "origprot", "era", "race_bw", "hisp")
CENS_COVARS <- c("age", "male", "edu", "apoe4", "cdrsb", "black", "hispanic", "era")
SCEN     <- c("S1 main outcome, main-trained", "S2 confirmed outcome, main-trained", "S3 confirmed outcome, confirmed-trained")
PAL      <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000")
# ---------------------------------------------------------------------------

tests <- data.frame(test = character(), result = character(), detail = character(), stringsAsFactors = FALSE)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s %s\n", if (ok) "PASS" else "FAIL", name, detail))
  tests <<- rbind(tests, data.frame(test = name, result = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE))
  invisible(ok)
}
rd <- function(x, k = 3) { nn <- vapply(x, is.numeric, logical(1)); x[nn] <- lapply(x[nn], round, k); x }

# Kaplan-Meier risk by h without survfit (fast enough for bootstrap loops)
km_fast <- function(time, ev, h) {
  o <- order(time); t <- time[o]; e <- ev[o]
  te <- t[e == 1 & t <= h]
  if (!length(te)) return(0)
  r <- rle(te)
  nr <- length(t) - findInterval(r$values, t, left.open = TRUE)
  1 - prod(1 - r$lengths / nr)
}

oc    <- readRDS(file.path(CFG$dir_out, "interim", "outcomes.rds"))
vis   <- readRDS(file.path(CFG$dir_out, "interim", "visits.rds"))
folds <- readRDS(file.path(CFG$dir_out, "interim", "folds.rds"))
o07   <- readRDS(file.path(CFG$dir_out, "interim", "oof_predictions.rds"))

d0 <- NULL; PM <- list(); PC <- list()
ev_df <- ov_df <- delta_df <- sub_df <- con_df <- NULL
bounds <- NULL; ipcw_tab <- NULL

section("0. Data, confirmed outcome and out-of-fold predictions", {
  d0 <- make_d0(oc, vis, folds, "in_set")
  d0$time_c <- (ifelse(!is.na(d0$t_dem_conf), as.numeric(d0$t_dem_conf), as.numeric(d0$last_obs)) - as.numeric(d0$t0)) / DAYS
  d0$ev_c   <- as.integer(!is.na(d0$t_dem_conf))
  cat("Modelling rows:", nrow(d0), "| rows with a main event:", sum(d0$ev), "| rows with a confirmed event:", sum(d0$ev_c), "\n")
  for (l in LM_ALL) {
    dl <- d0[d0$landmark == l, ]
    PM[[l]] <- run_cv(dl, FB)$P
    dl2 <- dl; dl2$time <- dl$time_c; dl2$ev <- dl$ev_c
    PC[[l]] <- run_cv(dl2, FB)$P
  }
  cat("Refitted Cox-B (Tier B) on the main and on the confirmed outcome at", paste(LM_ALL, collapse = ", "), "\n")
})

section("A1. Events: main versus confirmed dementia", {
  rows <- list()
  for (l in LM_ALL) for (h in HORIZONS) {
    dl <- d0[d0$landmark == l, ]
    m <- sum(dl$ev == 1 & dl$time <= h); cf <- sum(dl$ev_c == 1 & dl$time_c <= h)
    rows[[length(rows) + 1]] <- data.frame(landmark = l, horizon = h, n = nrow(dl), events_main = m, events_confirmed = cf,
                                           not_confirmed = m - cf, pct_not_confirmed = round(100 * (m - cf) / max(m, 1), 1))
  }
  ev_df <- do.call(rbind, rows)
  save_csv(ev_df, "13_events_main_vs_confirmed.csv"); show(ev_df, 20)
})

# ---- A2 and A3: metrics and subgroup O/E under the three scenarios --------------
section("A2. Metrics and subgroup O/E under the three scenarios (this is the slow step)", {
  set.seed(CFG$seed + 1500)
  ov_rows <- list(); dl_rows <- list(); sb_rows <- list(); cn_rows <- list()
  nH <- length(HORIZONS); nHS <- length(HS_SUB)
  for (l in LM_ALL) {
    dl <- d0[d0$landmark == l, ]; n <- nrow(dl)
    SC <- list(list(time = dl$time,   ev = dl$ev,   P = PM[[l]]),
               list(time = dl$time_c, ev = dl$ev_c, P = PM[[l]]),
               list(time = dl$time_c, ev = dl$ev_c, P = PC[[l]]))
    specs <- list()
    if (l %in% LM_SUB) for (v in SUBV) {
      lv <- as.character(dl[[v]])
      for (u in setdiff(unique(lv), c(NA, "Unknown", "missing"))) {
        inG <- !is.na(lv) & lv == u
        if (sum(inG) >= MIN_N) specs[[length(specs) + 1]] <- list(var = v, level = u, inG = inG)
      }
    }
    nsp <- length(specs)
    # point estimates
    pt <- array(NA_real_, c(3, nH, 3), dimnames = list(NULL, NULL, c("AUC", "sBrier", "OE")))
    ptev <- matrix(NA_integer_, 3, nH)
    for (s in 1:3) for (j in seq_len(nH)) {
      pt[s, j, ] <- met3(SC[[s]]$time, SC[[s]]$ev, SC[[s]]$P[, j], HORIZONS[j])
      ptev[s, j] <- sum(SC[[s]]$ev == 1 & SC[[s]]$time <= HORIZONS[j])
    }
    psub <- array(NA_real_, c(max(nsp, 1), 3, nHS)); psev <- array(NA_integer_, c(max(nsp, 1), 3, nHS))
    for (sp in seq_len(nsp)) for (s in 1:3) for (hh in seq_len(nHS)) {
      ig <- which(specs[[sp]]$inG); h <- HS_SUB[hh]; jh <- match(h, HORIZONS)
      psub[sp, s, hh] <- km_fast(SC[[s]]$time[ig], SC[[s]]$ev[ig], h) / mean(SC[[s]]$P[ig, jh], na.rm = TRUE)
      psev[sp, s, hh] <- sum(SC[[s]]$ev[ig] == 1 & SC[[s]]$time[ig] <= h)
    }
    # bootstrap: one resample per replicate, shared by all scenarios
    auc_b <- oe_b <- array(NA_real_, c(3, nH, B))
    sub_b <- array(NA_real_, c(max(nsp, 1), 3, nHS, B))
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      for (s in 1:3) for (j in seq_len(nH)) {
        m <- met3(SC[[s]]$time[idx], SC[[s]]$ev[idx], SC[[s]]$P[idx, j], HORIZONS[j])
        auc_b[s, j, b] <- m[["AUC"]]; oe_b[s, j, b] <- m[["OE"]]
      }
      for (sp in seq_len(nsp)) {
        ig <- idx[specs[[sp]]$inG[idx]]
        if (length(ig) < 2) next
        for (s in 1:3) for (hh in seq_len(nHS)) {
          jh <- match(HS_SUB[hh], HORIZONS)
          sub_b[sp, s, hh, b] <- km_fast(SC[[s]]$time[ig], SC[[s]]$ev[ig], HS_SUB[hh]) / mean(SC[[s]]$P[ig, jh], na.rm = TRUE)
        }
      }
    }
    for (s in 1:3) for (j in seq_len(nH)) {
      a <- q95(auc_b[s, j, ]); o <- q95(oe_b[s, j, ])
      ov_rows[[length(ov_rows) + 1]] <- data.frame(landmark = l, horizon = HORIZONS[j], scenario = SCEN[s], n = n,
        events_by_h = ptev[s, j], AUC = pt[s, j, "AUC"], AUC_lcl = a[1], AUC_ucl = a[2],
        sBrier = pt[s, j, "sBrier"], O_over_E = pt[s, j, "OE"], OE_lcl = o[1], OE_ucl = o[2], stringsAsFactors = FALSE)
    }
    for (j in seq_len(nH)) {
      d21 <- auc_b[2, j, ] - auc_b[1, j, ]; d32 <- auc_b[3, j, ] - auc_b[2, j, ]
      for (cmp in 1:2) {
        dd <- if (cmp == 1) d21 else d32; est <- if (cmp == 1) pt[2, j, "AUC"] - pt[1, j, "AUC"] else pt[3, j, "AUC"] - pt[2, j, "AUC"]
        qq <- q95(dd)
        dl_rows[[length(dl_rows) + 1]] <- data.frame(landmark = l, horizon = HORIZONS[j],
          comparison = if (cmp == 1) "AUC: confirmed outcome minus main outcome (same model)" else "AUC: confirmed-trained minus main-trained (both on confirmed outcome)",
          delta_AUC = est, lcl = qq[1], ucl = qq[2], stringsAsFactors = FALSE)
      }
    }
    for (sp in seq_len(nsp)) for (s in 1:3) for (hh in seq_len(nHS)) {
      ok <- psev[sp, s, hh] >= MIN_EV
      qq <- if (ok) q95(sub_b[sp, s, hh, ]) else c(NA_real_, NA_real_)
      sb_rows[[length(sb_rows) + 1]] <- data.frame(landmark = l, horizon = HS_SUB[hh], scenario = SCEN[s],
        variable = specs[[sp]]$var, level = specs[[sp]]$level, n = sum(specs[[sp]]$inG), events_by_h = psev[sp, s, hh],
        O_over_E = psub[sp, s, hh], lcl = qq[1], ucl = qq[2], stringsAsFactors = FALSE)
    }
    if (nsp) {
      ix <- function(v, u) { k <- which(vapply(specs, function(z) z$var == v && z$level == u, logical(1))); if (length(k) == 1) k else NA_integer_ }
      for (ct in list(c("origprot", "ADNI3", "ADNI1"), c("race_bw", "Black", "White"))) {
        a <- ix(ct[1], ct[2]); r <- ix(ct[1], ct[3]); if (is.na(a) || is.na(r)) next
        for (s in 1:3) for (hh in seq_len(nHS)) {
          ok <- psev[a, s, hh] >= MIN_EV && psev[r, s, hh] >= MIN_EV
          qq <- if (ok) q95(sub_b[a, s, hh, ] / sub_b[r, s, hh, ]) else c(NA_real_, NA_real_)
          cn_rows[[length(cn_rows) + 1]] <- data.frame(landmark = l, horizon = HS_SUB[hh], scenario = SCEN[s],
            contrast = paste(ct[2], "vs", ct[3]), ratio_of_OE = psub[a, s, hh] / psub[r, s, hh], lcl = qq[1], ucl = qq[2], stringsAsFactors = FALSE)
        }
      }
    }
    cat(sprintf("  %-4s done (%d rows, %d subgroup levels)\n", l, n, nsp))
  }
  ov_df <- do.call(rbind, ov_rows); delta_df <- do.call(rbind, dl_rows)
  sub_df <- do.call(rbind, sb_rows); con_df <- do.call(rbind, cn_rows)
  save_csv(rd(ov_df), "13_metrics_by_outcome.csv"); save_csv(rd(delta_df), "13_auc_differences.csv")
  save_csv(rd(sub_df), "13_subgroup_oe_by_outcome.csv"); save_csv(rd(con_df), "13_contrasts_by_outcome.csv")
})

section("A3. Results: whole cohort", {
  cat("AUC and O/E at 2 years, by landmark and scenario:\n")
  show(rd(ov_df[ov_df$horizon == 2, c("landmark", "scenario", "events_by_h", "AUC", "AUC_lcl", "AUC_ucl", "O_over_E", "OE_lcl", "OE_ucl")]), 20)
  cat("\nPaired AUC differences, 2 years:\n")
  show(rd(delta_df[delta_df$horizon == 2, ]), 20)
})

section("A4. Results: protocol era, race and ethnicity at baseline", {
  cat("Subgroup O/E, baseline, 2 years (intervals only where at least", MIN_EV, "events):\n")
  show(rd(sub_df[sub_df$landmark == "bl" & sub_df$horizon == 2, c("scenario", "variable", "level", "n", "events_by_h", "O_over_E", "lcl", "ucl")]), 60)
  cat("\nRatios of O/E (ADNI3 vs ADNI1, Black vs White):\n")
  show(rd(con_df[con_df$landmark %in% LM_SUB, ]), 30)
})

# ---- B. censoring bounds for Cox-B -----------------------------------------------
aud <- add_groups(oc[oc$in_audit, ])
PB <- list()
section("B0. Cox-B predictions from 07, aligned to the modelling rows", {
  for (l in LM_B) PB[[l]] <- o07$P[[l]][["Cox-B"]]
  cat("Using saved Cox-B out-of-fold predictions for", paste(LM_B, collapse = ", "), "\n")
})

section("B1. Observed-risk bounds and tipping fraction (audit set, observable subjects)", {
  rows <- list()
  for (l in LM_B) for (h in HS_B) {
    a  <- aud[aud$landmark == l, ]; cl <- a[[paste0("class_h", h)]]
    dl <- d0[d0$landmark == l, ];  dcl <- dl[[paste0("class_h", h)]]
    j  <- match(h, HORIZONS); P <- PB[[l]][, j]
    for (v in GRP_B) {
      lvA <- if (v == "all") rep("all", nrow(a)) else as.character(a[[v]])
      lvD <- if (v == "all") rep("all", nrow(dl)) else as.character(dl[[v]])
      for (u in setdiff(unique(lvA), c(NA, "Unknown", "missing"))) {
        ia <- which(lvA == u & cl != "administrative")
        if (length(ia) < MIN_N) next
        dem <- sum(cl[ia] == "dementia"); lost <- sum(cl[ia] == "lost")
        lnv <- sum(cl[ia] == "lost" & !a$has_followup[ia])
        idd <- which(lvD == u); io <- which(lvD == u & dcl != "administrative" & !is.na(P))
        pred <- mean(P[io])
        rows[[length(rows) + 1]] <- data.frame(landmark = l, horizon = h, variable = v, level = u,
          n_observable = length(ia), dementia = dem, lost = lost, lost_no_visit_after_landmark = lnv,
          risk_lower = dem / length(ia), risk_upper = (dem + lost) / length(ia),
          km_risk_modelling_set = km_fast(dl$time[idd], dl$ev[idd], h), mean_pred = pred,
          tipping_fraction = if (lost > 0) (pred * length(ia) - dem) / lost else NA_real_, stringsAsFactors = FALSE)
      }
    }
  }
  bounds <- do.call(rbind, rows)
  bounds$read <- ifelse(is.na(bounds$tipping_fraction), "no lost subjects",
                  ifelse(bounds$tipping_fraction < 0, "observed already above predicted",
                  ifelse(bounds$tipping_fraction > 1, "no censoring scenario reconciles", "censoring could reconcile")))
  save_csv(rd(bounds), "13_bounds_tipping_coxb.csv")
  cat("Cox-B, 2 years (tipping_fraction = share of lost subjects who must have converted for observed = predicted;\n",
      "the lost subjects are compared with the fraction among those followed in the km_risk column):\n", sep = "")
  show(rd(bounds[bounds$horizon == 2 & bounds$landmark == "bl", c("variable", "level", "n_observable", "dementia", "lost",
        "lost_no_visit_after_landmark", "risk_lower", "risk_upper", "km_risk_modelling_set", "mean_pred", "tipping_fraction", "read")]), 40)
})

section("B2. Covariate-adjusted IPCW observed risk versus Kaplan-Meier (Cox-B predictions)", {
  rows <- list()
  for (l in LM_B) {
    dl <- d0[d0$landmark == l, ]
    dl$black <- as.integer(dl$race == "Black"); dl$hispanic <- as.integer(dl$hisp == "Hispanic")
    for (f in c("age", "male", "edu", "apoe4", "cdrsb")) dl[[f]][is.na(dl[[f]])] <- median(dl[[f]], na.rm = TRUE)
    fit_c <- survival::coxph(stats::as.formula(paste("survival::Surv(time, 1 - ev) ~", paste(CENS_COVARS, collapse = " + "))), data = dl)
    risk <- predict(fit_c, type = "risk")
    bh <- survival::basehaz(fit_c, centered = TRUE); H0 <- stepfun(bh$time, c(0, bh$hazard))
    GT <- pmax(exp(-H0(dl$time) * risk), 0.05)
    for (h in HS_B) {
      Gh <- pmax(exp(-H0(h) * risk), 0.05); j <- match(h, HORIZONS)
      for (v in GRP_B) {
        lv <- if (v == "all") rep("all", nrow(dl)) else as.character(dl[[v]])
        for (u in setdiff(unique(lv), c(NA, "Unknown", "missing"))) {
          s <- which(lv == u)
          if (length(s) < MIN_N) next
          kmr <- km_fast(dl$time[s], dl$ev[s], h)
          ipr <- ipcw_risk(dl$time, dl$ev, GT, Gh, h, idx = lv == u & !is.na(lv))
          pm <- mean(PB[[l]][s, j], na.rm = TRUE)
          rows[[length(rows) + 1]] <- data.frame(landmark = l, horizon = h, variable = v, level = u, n = length(s),
            events_by_h = sum(dl$ev[s] == 1 & dl$time[s] <= h), km_risk = kmr, ipcw_risk = ipr, mean_pred = pm,
            OE_km = kmr / pm, OE_ipcw = ipr / pm, stringsAsFactors = FALSE)
        }
      }
    }
  }
  ipcw_tab <- do.call(rbind, rows)
  save_csv(rd(ipcw_tab), "13_ipcw_vs_km_coxb.csv")
  show(rd(ipcw_tab[ipcw_tab$landmark == "bl" & ipcw_tab$horizon == 2, c("variable", "level", "n", "events_by_h", "km_risk", "ipcw_risk", "mean_pred", "OE_km", "OE_ipcw")]), 40)
  cat("\nIf OE_ipcw stays close to OE_km, covariate-explained censoring does not drive the miscalibration.\n")
})

section("C. Figure: era calibration under the three outcome definitions", {
  fdir <- file.path(CFG$dir_out, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  ERAS <- c("ADNI1", "ADNI2", "ADNI3", "ADNI4")
  ci_v <- function(x, l, u, lo, hi, col) {
    if (!is.finite(l) || !is.finite(u)) return(invisible())
    segments(x, max(l, lo), x, min(u, hi), col = col, lwd = 1.5)
    if (l < lo) points(x, lo, pch = "v", col = col, cex = 0.7)
    if (u > hi) points(x, hi, pch = "^", col = col, cex = 0.7)
  }
  fig <- function() {
    par(mfrow = c(1, 2), mar = c(5, 4.2, 3, 1))
    for (l in LM_SUB) {
      plot(NA, xlim = c(0.5, length(ERAS) + 0.5), ylim = c(0.2, 3), log = "y", xaxt = "n", yaxt = "n",
           xlab = "Protocol era", ylab = "Observed / expected (2-year horizon)", main = paste("Landmark", l))
      axis(1, at = seq_along(ERAS), labels = ERAS)
      axis(2, at = c(0.25, 0.5, 0.67, 1, 1.5, 2, 3), labels = c("0.25", "0.5", "0.67", "1", "1.5", "2", "3"), las = 1)
      abline(h = 1, lty = 2, col = "grey40")
      for (s in 1:3) for (k in seq_along(ERAS)) {
        r <- sub_df[sub_df$landmark == l & sub_df$horizon == 2 & sub_df$variable == "origprot" & sub_df$level == ERAS[k] & sub_df$scenario == SCEN[s], ]
        if (nrow(r) != 1 || r$events_by_h < MIN_EV) next
        x <- k + (s - 2) * 0.2
        ci_v(x, r$lcl, r$ucl, 0.2, 3, PAL[s]); points(x, r$O_over_E, pch = 16, col = PAL[s])
      }
      legend("topright", legend = c("Main outcome", "Confirmed outcome, main-trained", "Confirmed outcome, confirmed-trained"),
             col = PAL[1:3], pch = 16, bty = "n", cex = 0.65)
    }
  }
  ok <- tryCatch({
    grDevices::png(file.path(fdir, "13_era_oe_by_outcome.png"), width = 9, height = 4.6, units = "in", res = 200); fig(); grDevices::dev.off()
    grDevices::pdf(file.path(fdir, "13_era_oe_by_outcome.pdf"), width = 9, height = 4.6); fig(); grDevices::dev.off(); TRUE
  }, error = function(e) { while (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE); cat("  figure not written:", conditionMessage(e), "\n"); FALSE })
  assign("fig_ok13", ok, envir = .GlobalEnv)
  cat("Figure:", if (ok) "13_era_oe_by_outcome.png and .pdf" else "none", "\n")
})

section("D. Checks", {
  ok_pred <- !is.null(d0) && length(PM) == length(LM_ALL)
  check("data and predictions were built", ok_pred)
  if (ok_pred) {
    key <- function(x) paste(x$RID, x$landmark)
    check("modelling rows line up with the rows saved by 07", identical(key(d0), key(o07$d)))
    same <- all(vapply(LM_ALL, function(l) isTRUE(all.equal(unname(PM[[l]]), unname(o07$P[[l]][["Cox-B"]]), tolerance = 1e-8)), logical(1)))
    check("refitting Cox-B on the main outcome reproduces the 07 predictions", same)
    dl <- d0[d0$landmark == "bl", ]
    check("km_fast reproduces survfit Kaplan-Meier risk (within 1e-9)", abs(km_fast(dl$time, dl$ev, 2) - km_risk(dl$time, dl$ev, 2)) < 1e-9)
    check("whole-cohort baseline events at 2 years are 201", sum(dl$ev == 1 & dl$time <= 2) == 201, sprintf("(%d)", sum(dl$ev == 1 & dl$time <= 2)))
    check("confirmed events never exceed main events", all(ev_df$events_confirmed <= ev_df$events_main))
    e3 <- tryCatch(readr::read_csv(file.path(CFG$dir_out, "audit", "03_events_by_horizon.csv"), show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(e3)) {
      m <- merge(ev_df, e3[, c("landmark", "horizon", "dementia_confirmed")], by = c("landmark", "horizon"))
      check("confirmed event counts match 03", nrow(m) == nrow(ev_df) && all(m$events_confirmed == m$dementia_confirmed),
            sprintf("(%d of %d rows equal)", sum(m$events_confirmed == m$dementia_confirmed), nrow(ev_df)))
    } else check("03 events file found", FALSE)
    a1 <- ov_df[ov_df$landmark == "bl" & ov_df$horizon == 2 & ov_df$scenario == SCEN[1], ]
    check("S1 AUC at baseline, 2 years matches the 07 report (0.852)", nrow(a1) == 1 && abs(a1$AUC - 0.852) < 0.003, sprintf("(%.3f)", a1$AUC))
    r8 <- tryCatch(readr::read_csv(file.path(CFG$dir_out, "audit", "08_oe_over_time.csv"), show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(r8)) {
      x <- r8[r8$variable == "all" & r8$model == "Cox-B" & r8$landmark == "bl" & r8$horizon == 2, ]
      check("S1 whole-cohort O/E at baseline, 2 years matches 08", nrow(x) == 1 && abs(x$O_over_E - a1$O_over_E) < 0.002, sprintf("(%.3f vs %.3f)", a1$O_over_E, x$O_over_E))
      y <- sub_df[sub_df$landmark == "bl" & sub_df$horizon == 2 & sub_df$scenario == SCEN[1] & sub_df$variable == "origprot" & sub_df$level == "ADNI3", ]
      z <- r8[r8$variable == "origprot" & r8$level == "ADNI3" & r8$model == "Cox-B" & r8$landmark == "bl" & r8$horizon == 2, ]
      check("S1 ADNI3 O/E at baseline, 2 years matches 08", nrow(y) == 1 && nrow(z) == 1 && abs(y$O_over_E - z$O_over_E) < 0.002, sprintf("(%.3f vs %.3f)", y$O_over_E, z$O_over_E))
    } else check("08 results found", FALSE)
    ok_ci <- ov_df$AUC >= ov_df$AUC_lcl - 0.05 & ov_df$AUC <= ov_df$AUC_ucl + 0.05
    check("AUC point estimates sit inside their bootstrap intervals (within 0.05)", all(ok_ci, na.rm = TRUE), sprintf("(%d of %d)", sum(ok_ci, na.rm = TRUE), sum(!is.na(ok_ci))))
    check("risk bounds are ordered and inside 0 to 1", !is.null(bounds) && all(bounds$risk_lower <= bounds$risk_upper & bounds$risk_upper <= 1))
    bb <- bounds[bounds$variable == "all" & bounds$landmark == "bl" & bounds$horizon == 2, ]
    check("Kaplan-Meier risk lies between the bounds for the whole cohort (within 0.05)", nrow(bb) == 1 && bb$km_risk_modelling_set >= bb$risk_lower - 0.05 && bb$km_risk_modelling_set <= bb$risk_upper + 0.05)
    check("IPCW risks are between 0 and 1", !is.null(ipcw_tab) && all(ipcw_tab$ipcw_risk >= 0 & ipcw_tab$ipcw_risk <= 1, na.rm = TRUE))
    check("figure written", isTRUE(get0("fig_ok13", ifnotfound = FALSE)) && file.exists(file.path(CFG$dir_out, "figures", "13_era_oe_by_outcome.png")))
  }
  cat(sprintf("\n%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  save_csv(tests, "13_tests.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
