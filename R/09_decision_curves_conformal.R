# 09_decision_curves_conformal.R ---------------------------------------------
# Two things from the out-of-fold predictions saved by 07 (no models are refit):
#
# 1. Decision curves. Net benefit of acting on a model's predicted risk at each
#    risk threshold, against treating everyone and treating no one, with
#    censoring handled by Kaplan-Meier risk inside the treated group (Vickers
#    et al. 2008). A threshold pt stands for the odds of harm to benefit you
#    accept: pt = 0.20 means one unnecessary intervention per four true cases.
#    Bootstrap intervals are over subjects. Subgroup curves (point estimates
#    only) show what the era miscalibration found in 08 costs in net benefit.
#
# 2. Conformal prediction sets for "dementia within h years" (yes / no).
#    Cross-conformal: each fold is calibrated on the other four folds. Subjects
#    censored before h have unknown status and cannot be scored, so calibration
#    and coverage use inverse-probability-of-censoring weights (1/G(T) for cases,
#    1/G(h) for controls); unweighted ("naive") versions are shown beside them.
#    Both a marginal and a class-conditional (Mondrian) version are run, because
#    events are rare. Coverage is then tracked by subgroup.
#
# Reporting rules: subgroup levels with fewer than MIN_N subjects are not
# analysed; levels with fewer than MIN_EV events by the horizon are listed but
# not given intervals. Death is treated as censoring, as in 04 to 08.
#
# Writes: output/audit/09_*.csv, 09_decision_curves_conformal_report.txt
#         output/figures/09_*.png (aggregate results only)
# Run:    source("R/09_decision_curves_conformal.R")   (after 07; a few minutes)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
if (!requireNamespace("survival", quietly = TRUE))
  stop("Package 'survival' is missing. It ships with R; run install.packages('survival').")
K <- 5; HORIZONS <- CFG$horizons_years
source(if (file.exists("R/utils_eval.R")) "R/utils_eval.R" else "utils_eval.R")

report <- file.path(CFG$dir_out, "audit", "09_decision_curves_conformal_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Decision curves and conformal prediction | run:", format(Sys.time()), "\n")

# ---- settings -------------------------------------------------------------
B        <- 500                                   # bootstrap replicates, decision curves
B_CF     <- 300                                   # bootstrap replicates, conformal coverage
MIN_N    <- 30
MIN_EV   <- 5
ALPHA    <- 0.10                                  # target coverage 90%
THR      <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
VARS     <- c("sex", "race_bw", "hisp", "edu_grp", "apoe4_grp", "age_grp", "origprot")
DC_SUB_MODEL <- "Cox-B"
CF_MODELS    <- c("Cox-B", "GLMNET-C", "RSF-C")
FIG_H <- 2
# ---------------------------------------------------------------------------

tests <- data.frame(test = character(), result = character(), detail = character(), stringsAsFactors = FALSE)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s %s\n", if (ok) "PASS" else "FAIL", name, detail))
  tests <<- rbind(tests, data.frame(test = name, result = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE))
  invisible(ok)
}

# ---- helpers --------------------------------------------------------------
# Kaplan-Meier risk by h (fast, tie-aware; agrees with survfit, see checks).
km_fast <- function(time, ev, h) {
  n <- length(time); if (n == 0) return(NA_real_)
  idx <- which(ev == 1 & time <= h)
  if (!length(idx)) return(0)
  ts <- sort(time); ut <- sort(unique(time[idx]))
  d  <- tabulate(match(time[idx], ut), length(ut))
  nr <- n - findInterval(ut, ts, left.open = TRUE)          # number with time >= u
  1 - prod(1 - d / nr)
}
# Net benefit of treating everyone with predicted risk >= pt. Returns c(nb, share treated).
nb_fun <- function(time, ev, p, h, pt) {
  n <- length(time); g <- which(p >= pt)
  if (!length(g)) return(c(nb = 0, share = 0))
  r <- km_fast(time[g], ev[g], h); s <- length(g) / n
  c(nb = s * r - s * (1 - r) * pt / (1 - pt), share = s)
}
nb_all <- function(time, ev, h, pt) { r <- km_fast(time, ev, h); r - (1 - r) * pt / (1 - pt) }

# Weighted conformal quantile (test point carries weight 1 of mean-1 weights).
wq <- function(s, w, alpha) {
  if (!length(s)) return(Inf)
  o <- order(s); s <- s[o]; w <- w[o]; w <- w / mean(w)
  k <- which(cumsum(w) / (sum(w) + 1) >= 1 - alpha)[1]
  if (is.na(k)) Inf else s[k]
}
# Cross-conformal sets: y is 1 (case by h), 0 (control past h), NA (censored before h).
conformal_cv <- function(p, y, w, fold, alpha, classcond, weighted) {
  n <- length(p); in1 <- in0 <- logical(n)
  ww <- if (weighted) w else rep(1, n)
  for (f in sort(unique(fold))) {
    te <- fold == f; ca <- !te & !is.na(y)
    if (classcond) {
      c1 <- ca & y == 1; c0 <- ca & y == 0
      q1 <- wq(1 - p[c1], ww[c1], alpha); q0 <- wq(p[c0], ww[c0], alpha)
    } else {
      q1 <- q0 <- wq(ifelse(y[ca] == 1, 1 - p[ca], p[ca]), ww[ca], alpha)
    }
    in1[te] <- (1 - p[te]) <= q1; in0[te] <- p[te] <= q0
  }
  list(in1 = in1, in0 = in0)
}
wmean <- function(x, w) if (length(x) == 0 || sum(w) == 0) NA_real_ else sum(w * x) / sum(w)
rd <- function(x) { nn <- vapply(x, is.numeric, logical(1)); x[nn] <- lapply(x[nn], round, 3); x }

dat <- readRDS(file.path(CFG$dir_out, "interim", "oof_predictions.rds"))
D <- dat$d; Pm <- dat$P
stopifnot(identical(as.numeric(dat$horizons), as.numeric(HORIZONS)))
LM <- CFG$landmarks

# subgroup level specifications for one landmark's rows
make_specs <- function(dl) {
  n <- nrow(dl); specs <- list(list(var = "all", level = "all", idx = seq_len(n)))
  for (v in VARS) {
    lv <- as.character(dl[[v]])
    for (u in setdiff(unique(lv), c(NA, "Unknown", "missing")))
      if (sum(lv == u, na.rm = TRUE) >= MIN_N) specs[[length(specs) + 1]] <- list(var = v, level = u, idx = which(lv == u))
  }
  specs
}

# ---- 1. decision curves ---------------------------------------------------
dc_list <- list(); dcs_list <- list()
section("1. Decision curves: whole cohort, every landmark, horizon and model", {
  for (l in LM) {
    dl <- D[D$landmark == l, ]; n <- nrow(dl); Pl <- Pm[[l]]; mods <- names(Pl)
    set.seed(CFG$seed + 900 + match(l, LM))
    bidx <- replicate(B, sample.int(n, n, replace = TRUE), simplify = FALSE)
    for (j in seq_along(HORIZONS)) {
      h <- HORIZONS[j]; evn <- sum(dl$ev == 1 & dl$time <= h)
      nba   <- sapply(THR, function(pt) nb_all(dl$time, dl$ev, h, pt))
      nba_b <- t(sapply(bidx, function(ix) sapply(THR, function(pt) nb_all(dl$time[ix], dl$ev[ix], h, pt))))
      for (m in mods) {
        p <- Pl[[m]][, j]
        est <- sapply(THR, function(pt) nb_fun(dl$time, dl$ev, p, h, pt))          # 2 x T
        nb_b <- matrix(NA_real_, B, length(THR))
        for (b in seq_len(B)) { ix <- bidx[[b]]; pb <- p[ix]
          for (k in seq_along(THR)) nb_b[b, k] <- nb_fun(dl$time[ix], dl$ev[ix], pb, h, THR[k])[1] }
        dlt <- est[1, ] - pmax(nba, 0); dlt_b <- nb_b - pmax(nba_b, 0)
        ci <- t(apply(dlt_b, 2, q95))
        dc_list[[length(dc_list) + 1]] <- data.frame(landmark = l, model = m, horizon = h, threshold = THR,
          n = n, events_by_h = evn, NB_model = est[1, ], NB_treat_all = nba, NB_best_default = pmax(nba, 0),
          delta_NB = dlt, lcl = ci[, 1], ucl = ci[, 2], share_treated = est[2, ], stringsAsFactors = FALSE)
      }
    }
    cat("  ", l, ": done\n")
  }
  dc_df <- rd(do.call(rbind, dc_list)); save_csv(dc_df, "09_decision_curve.csv")
  assign("dc_df", dc_df, envir = .GlobalEnv)
  cat("\nCox-B at baseline, 2-year horizon (delta_NB = net benefit minus the better of treat-all and treat-none):\n")
  show(dc_df[dc_df$landmark == "bl" & dc_df$model == "Cox-B" & dc_df$horizon == 2,
             c("threshold", "NB_model", "NB_treat_all", "delta_NB", "lcl", "ucl", "share_treated")], 20)
  # where does the model beat the default with an interval above 0?
  win <- dc_df[!is.na(dc_df$lcl) & dc_df$lcl > 0, ]
  cat("\nThresholds at which the interval for delta_NB lies above 0 (model clearly beats the default), by landmark, horizon, model:\n")
  if (nrow(win)) {
    agg <- aggregate(threshold ~ landmark + horizon + model, win, function(x) paste0(min(x), "-", max(x), " (", length(x), " of ", length(THR), ")"))
    show(agg[order(agg$landmark, agg$horizon, agg$model), ], 80)
  } else cat("  none\n")
})

section("2. Decision curves by subgroup", {
  sub_model <- DC_SUB_MODEL
  for (l in LM) {
    dl <- D[D$landmark == l, ]; Pl <- Pm[[l]]; specs <- make_specs(dl)
    for (j in seq_along(HORIZONS)) {
      h <- HORIZONS[j]; p <- Pl[[sub_model]][, j]
      for (s in specs) {
        i <- s$idx; evn <- sum(dl$ev[i] == 1 & dl$time[i] <= h)
        if (evn < MIN_EV) next
        nba <- sapply(THR, function(pt) nb_all(dl$time[i], dl$ev[i], h, pt))
        est <- sapply(THR, function(pt) nb_fun(dl$time[i], dl$ev[i], p[i], h, pt))
        dcs_list[[length(dcs_list) + 1]] <- data.frame(landmark = l, model = sub_model, horizon = h, variable = s$var, level = s$level,
          n = length(i), events_by_h = evn, threshold = THR, NB_model = est[1, ], NB_treat_all = nba,
          delta_NB = est[1, ] - pmax(nba, 0), share_treated = est[2, ], stringsAsFactors = FALSE)
      }
    }
  }
  dcs_df <- rd(do.call(rbind, dcs_list)); save_csv(dcs_df, "09_decision_curve_by_subgroup.csv")
  assign("dcs_df", dcs_df, envir = .GlobalEnv)
  cat("Point estimates only (no intervals); levels with fewer than", MIN_EV, "events by the horizon are left out.\n")
  cat("\nBaseline, 2-year horizon, threshold 0.20, by protocol era, sex and APOE4:\n")
  x <- dcs_df[dcs_df$landmark == "bl" & dcs_df$horizon == 2 & dcs_df$threshold == 0.20 &
                dcs_df$variable %in% c("all", "origprot", "sex", "apoe4_grp"),
              c("variable", "level", "n", "events_by_h", "NB_model", "NB_treat_all", "delta_NB", "share_treated")]
  show(x, 30)
})

# ---- 3. conformal prediction sets ---------------------------------------------
cf_over <- list(); cf_sub <- list()
METHODS <- data.frame(method = c("marginal-IPCW", "class-cond-IPCW", "marginal-naive", "class-cond-naive"),
                      classcond = c(FALSE, TRUE, FALSE, TRUE), weighted = c(TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)
section("3. Conformal prediction sets (target coverage 90%)", {
  set.seed(CFG$seed + 950)
  for (l in LM) {
    dl <- D[D$landmark == l, ]; n <- nrow(dl); Pl <- Pm[[l]]; specs <- make_specs(dl)
    Gf <- cens_G(dl$time, dl$ev)
    for (j in seq_along(HORIZONS)) {
      h <- HORIZONS[j]; Gh <- Gf(h)
      y <- ifelse(dl$ev == 1 & dl$time <= h, 1, ifelse(dl$time > h, 0, NA))
      w <- ifelse(is.na(y), NA, ifelse(y == 1, 1 / Gf(dl$time), 1 / Gh))
      kn <- !is.na(y)
      for (m in intersect(CF_MODELS, names(Pl))) {
        p <- Pl[[m]][, j]
        for (mi in seq_len(nrow(METHODS))) {
          cs <- conformal_cv(p, y, w, dl$fold, ALPHA, METHODS$classcond[mi], METHODS$weighted[mi])
          cov <- ifelse(y == 1, cs$in1, cs$in0)
          sz1 <- cs$in1 != cs$in0
          cf_over[[length(cf_over) + 1]] <- data.frame(landmark = l, model = m, horizon = h, method = METHODS$method[mi],
            n_known = sum(kn), n_cases = sum(y == 1, na.rm = TRUE), n_unknown = sum(!kn),
            cov_all_ipcw = wmean(cov[kn], w[kn]), cov_events_ipcw = wmean(cov[which(y == 1)], w[which(y == 1)]),
            cov_nonevents_ipcw = wmean(cov[which(y == 0)], w[which(y == 0)]),
            cov_all_unw = mean(cov[kn]), cov_events_unw = mean(cov[which(y == 1)]), cov_nonevents_unw = mean(cov[which(y == 0)]),
            frac_singleton = mean(sz1), frac_both = mean(cs$in1 & cs$in0), frac_empty = mean(!cs$in1 & !cs$in0),
            stringsAsFactors = FALSE)
          # subgroup coverage for the two IPCW methods
          if (METHODS$weighted[mi]) for (s in specs) {
            i <- s$idx[kn[s$idx]]; nc <- sum(y[i] == 1)
            pt <- wmean(cov[i], w[i])
            ci <- c(NA_real_, NA_real_)
            if (nc >= MIN_EV && length(i) >= MIN_N) {
              bb <- replicate(B_CF, { k <- sample(i, length(i), replace = TRUE); wmean(cov[k], w[k]) })
              ci <- q95(bb)
            }
            cf_sub[[length(cf_sub) + 1]] <- data.frame(landmark = l, model = m, horizon = h, method = METHODS$method[mi],
              variable = s$var, level = s$level, n_subjects = length(s$idx), n_known = length(i), n_cases = nc,
              coverage_ipcw = pt, lcl = ci[1], ucl = ci[2], cov_events_ipcw = wmean(cov[i[y[i] == 1]], w[i[y[i] == 1]]),
              reported = nc >= MIN_EV, frac_singleton = mean(sz1[s$idx]), stringsAsFactors = FALSE)
          }
        }
      }
    }
    cat("  ", l, ": done\n")
  }
  cfo <- rd(do.call(rbind, cf_over)); cfs <- rd(do.call(rbind, cf_sub))
  save_csv(cfo, "09_conformal_overall.csv"); save_csv(cfs, "09_conformal_by_subgroup.csv")
  assign("cfo", cfo, envir = .GlobalEnv); assign("cfs", cfs, envir = .GlobalEnv)
  cat("\nWhole cohort, Cox-B, 2-year horizon (coverage is IPCW-weighted unless marked unw; target 0.90):\n")
  show(cfo[cfo$model == "Cox-B" & cfo$horizon == 2,
           c("landmark", "method", "n_known", "n_cases", "n_unknown", "cov_all_ipcw", "cov_events_ipcw", "cov_nonevents_ipcw",
             "cov_events_unw", "cov_nonevents_unw", "frac_singleton", "frac_both", "frac_empty")], 60)
})

section("4. Conformal coverage by subgroup (class-conditional, IPCW; Cox-B, 2-year horizon)", {
  x <- cfs[cfs$model == "Cox-B" & cfs$horizon == 2 & cfs$method == "class-cond-IPCW" & cfs$landmark == "bl",
           c("variable", "level", "n_known", "n_cases", "coverage_ipcw", "lcl", "ucl", "cov_events_ipcw", "frac_singleton", "reported")]
  show(x, 60)
  r <- cfs[cfs$reported & !is.na(cfs$lcl) & cfs$method == "class-cond-IPCW", ]
  below <- r[r$ucl < 1 - ALPHA, ]; above <- r[r$lcl > 1 - ALPHA, ]
  cat(sprintf("\nAcross landmarks, horizons and models (class-conditional): %d of %d subgroup intervals lie wholly below 0.90, %d wholly above it.\n",
              nrow(below), nrow(r), nrow(above)))
  if (nrow(below)) { cat("Below 0.90:\n"); show(below[, c("landmark", "model", "horizon", "variable", "level", "n_cases", "coverage_ipcw", "lcl", "ucl")], 60) }
  cat("\nIntervals reflect subject resampling at the evaluation step only, not uncertainty in the calibration quantile.\n")
  cat("Marginal-method subgroup rows are in 09_conformal_by_subgroup.csv.\n")
})

# ---- 5. figures -----------------------------------------------------------------
section("5. Figures", {
  fdir <- file.path(CFG$dir_out, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  ok1 <- tryCatch({
    f <- file.path(fdir, "09_decision_curves.png")
    grDevices::png(f, width = 2000, height = 600, res = 150)
    par(mfrow = c(1, 4), mar = c(4, 4, 3, 1))
    mods <- c("Cox-B", "GLMNET-C", "RSF-C"); cols <- grDevices::hcl.colors(3, "Dark 3")
    for (l in LM) {
      d <- dc_df[dc_df$landmark == l & dc_df$horizon == FIG_H, ]
      ymax <- max(d$NB_model, d$NB_treat_all, na.rm = TRUE)
      plot(NA, xlim = range(THR), ylim = c(-0.02, ymax * 1.05), xlab = "Risk threshold", ylab = "Net benefit",
           main = paste0(l, " (", FIG_H, "-year)"), cex.main = 0.9)
      abline(h = 0, col = "grey40")
      x <- d[d$model == "Cox-B", ]; lines(x$threshold, x$NB_treat_all, col = "grey55", lty = 2, lwd = 1.5)
      for (i in seq_along(mods)) { x <- d[d$model == mods[i], ]; if (nrow(x)) { lines(x$threshold, x$NB_model, col = cols[i], lwd = 1.8); points(x$threshold, x$NB_model, col = cols[i], pch = 19, cex = 0.5) } }
      if (l == LM[1]) legend("topright", legend = c("Treat all", mods), col = c("grey55", cols), lty = c(2, 1, 1, 1), bty = "n", cex = 0.7)
    }
    grDevices::dev.off(); TRUE
  }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); cat("  figure 1 not written:", conditionMessage(e), "\n"); FALSE })
  ok2 <- tryCatch({
    f <- file.path(fdir, "09_conformal_coverage_by_subgroup.png")
    grDevices::png(f, width = 1800, height = 1100, res = 150)
    par(mfrow = c(1, 2), mar = c(4, 9, 3, 1))
    for (meth in c("marginal-IPCW", "class-cond-IPCW")) {
      d <- cfs[cfs$model == "Cox-B" & cfs$horizon == FIG_H & cfs$landmark == "bl" & cfs$method == meth & cfs$reported, ]
      d <- d[nrow(d):1, ]; yy <- seq_len(nrow(d)); lab <- paste0(d$variable, ": ", d$level, " (", d$n_cases, " ev)")
      plot(NA, xlim = c(0.5, 1), ylim = c(0.5, nrow(d) + 0.5), yaxt = "n", xlab = "Coverage (IPCW-weighted)", ylab = "",
           main = paste0(meth, ", baseline, ", FIG_H, "-year"), cex.main = 0.9)
      axis(2, at = yy, labels = lab, las = 1, cex.axis = 0.6); abline(v = 1 - ALPHA, lty = 2, col = "grey40")
      segments(pmax(d$lcl, 0.5), yy, pmin(d$ucl, 1), yy, col = "#3b6ea5", lwd = 1.5)
      points(d$coverage_ipcw, yy, pch = 19, col = "#3b6ea5", cex = 0.7)
    }
    grDevices::dev.off(); TRUE
  }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); cat("  figure 2 not written:", conditionMessage(e), "\n"); FALSE })
  cat("Figures written to output/figures/ :", c("09_decision_curves.png", "09_conformal_coverage_by_subgroup.png")[c(ok1, ok2)], "\n")
})

# ---- 6. checks ------------------------------------------------------------------
section("6. Checks", {
  # fast KM agrees with survfit
  set.seed(1); dd <- D[D$landmark == "bl", ]; ii <- sample(nrow(dd), 400)
  kf <- sapply(HORIZONS, function(h) km_fast(dd$time[ii], dd$ev[ii], h))
  ks <- 1 - summary(survival::survfit(survival::Surv(dd$time[ii], dd$ev[ii]) ~ 1), times = HORIZONS, extend = TRUE)$surv
  check("fast Kaplan-Meier agrees with survfit at every horizon", all(abs(kf - ks) < 1e-9), sprintf("(largest difference %.2e)", max(abs(kf - ks))))
  # regression against 08
  r8 <- tryCatch(readr::read_csv(file.path(CFG$dir_out, "audit", "08_oe_over_time.csv"), show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(r8)) {
    a <- r8[r8$variable == "all" & r8$model == "Cox-B", ]
    ok <- sapply(seq_len(nrow(a)), function(k) {
      d1 <- D[D$landmark == a$landmark[k], ]; abs(km_fast(d1$time, d1$ev, a$horizon[k]) - a$km_observed[k]) < 0.002 })
    check("whole-cohort observed risk reproduces 08 for every landmark and horizon (within 0.002)", all(ok), sprintf("(%d rows)", nrow(a)))
  } else check("08 results found for the regression check", FALSE)
  # perfect predictor: net benefit equals prevalence at any threshold
  set.seed(2); ev <- rbinom(2000, 1, 0.3); tm <- ifelse(ev == 1, runif(2000, 0.1, 1.9), 5)
  nbp <- nb_fun(tm, ev, as.numeric(ev), 2, 0.2)[1]
  check("synthetic: a perfect predictor has net benefit equal to prevalence", abs(nbp - mean(ev)) < 1e-9, sprintf("(NB %.4f, prevalence %.4f)", nbp, mean(ev)))
  # treating everyone reproduces treat-all
  check("synthetic: threshold 0 treats everyone, so net benefit equals treat-all",
        abs(nb_fun(tm, ev, rep(0.5, 2000), 2, 0.2)[1] - nb_all(tm, ev, 2, 0.2)) < 1e-9)
  # negative control: shuffled outcomes give no net benefit over the default
  dd <- D[D$landmark == "bl", ]; set.seed(3); sh <- sample(nrow(dd)); p0 <- Pm[["bl"]][["Cox-B"]][, 2]
  dl0 <- sapply(THR, function(pt) nb_fun(dd$time[sh], dd$ev[sh], p0, 2, pt)[1] - max(nb_all(dd$time[sh], dd$ev[sh], 2, pt), 0))
  check("negative control: with shuffled outcomes the model never beats the default by more than 0.02", max(dl0) < 0.02, sprintf("(largest gain %.4f)", max(dl0)))
  # conformal synthetic: exchangeable, no censoring, true probabilities
  set.seed(4); n <- 6000; pp <- rbeta(n, 0.8, 2.5); yy <- rbinom(n, 1, pp); ff <- sample(rep(1:5, length.out = n))
  cm <- conformal_cv(pp, yy, rep(1, n), ff, ALPHA, FALSE, FALSE); cc <- conformal_cv(pp, yy, rep(1, n), ff, ALPHA, TRUE, FALSE)
  cvm <- mean(ifelse(yy == 1, cm$in1, cm$in0)); cv1 <- mean(cc$in1[yy == 1]); cv0 <- mean(cc$in0[yy == 0])
  check("synthetic: marginal conformal covers at least 88% and class-conditional covers each class at least 88%",
        cvm > 0.88 && cv1 > 0.88 && cv0 > 0.88, sprintf("(marginal %.3f, events %.3f, non-events %.3f)", cvm, cv1, cv0))
  # real data: standard conformal on the known-status subjects holds per class
  x <- cfo[cfo$landmark == "bl" & cfo$model == "Cox-B" & cfo$horizon == 2 & cfo$method == "class-cond-naive", ]
  check("real data: class-conditional naive coverage is at least 0.85 for events and non-events (bl, Cox-B, 2 years)",
        nrow(x) == 1 && x$cov_events_unw >= 0.85 && x$cov_nonevents_unw >= 0.85, sprintf("(events %.3f, non-events %.3f)", x$cov_events_unw, x$cov_nonevents_unw))
  check("saved predictions line up with the outcome rows at every landmark",
        all(sapply(LM, function(l) nrow(D[D$landmark == l, ]) == nrow(Pm[[l]][[1]]))))
  cat(sprintf("\n%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  save_csv(tests, "09_tests.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
