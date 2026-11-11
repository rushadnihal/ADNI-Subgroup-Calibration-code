# 06_censoring_power.R -----------------------------------------------------
# Two questions about the subgroup calibration results from 04 and 05.
#
# A. Could loss to follow-up explain the miscalibration we see (ADNI3, Black participants)?
#    A1. Bounds: observed 2- and 3-year risk if none / all of the lost subjects had converted,
#        plus the TIPPING FRACTION: the share of lost subjects who would need to have converted
#        for observed risk to equal predicted risk. Uses the audit set (includes subjects with no
#        visit after the landmark). Assumption-free.
#    A2. Covariate-adjusted IPCW observed risk, with the censoring model including age, sex,
#        education, APOE4, CDR-SB, Black, Hispanic and protocol era. Modelling set only, so it
#        corrects for loss after a subject returns, not for immediate dropout. Assumes censoring
#        depends only on those covariates.
# B. Power: the smallest observed/expected ratio each subgroup could detect with 80% power
#    (calibration-in-the-large, normal approximation to the Poisson-binomial, complete follow-up only).
#
# Predictions are the out-of-fold Tier T0 and Tier B Cox predictions on the 04 folds.
# Death is treated as censoring. Counts and metrics only in the CSVs.
#
# Writes: output/audit/06_*.csv, 06_censoring_power_report.txt
# Run:    source("R/06_censoring_power.R")   (after 05; about 1 to 2 minutes)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
if (!requireNamespace("survival", quietly = TRUE))
  stop("Package 'survival' is missing. It ships with R; run install.packages('survival').")
K <- 5; HORIZONS <- CFG$horizons_years
source(if (file.exists("R/utils_eval.R")) "R/utils_eval.R" else "utils_eval.R")

report <- file.path(CFG$dir_out, "audit", "06_censoring_power_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Censoring sensitivity and power | run:", format(Sys.time()), "\n")

LMS      <- c("bl", "m12")
HS       <- HORIZONS[HORIZONS >= 2]
GRP_VARS <- c("era", "race_bw", "hisp")
MIN_N    <- 20
TIERS    <- list(T0 = c("age", "apoe4", "cdrsb"),
                 B  = c("age", "male", "edu", "apoe4", "cdrsb", "adas13", "ravlt_immediate", "faq"))
CENS_COVARS <- c("age", "male", "edu", "apoe4", "cdrsb", "black", "hispanic", "era")

tests <- data.frame(test = character(), result = character(), detail = character(), stringsAsFactors = FALSE)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s %s\n", if (ok) "PASS" else "FAIL", name, detail))
  tests <<- rbind(tests, data.frame(test = name, result = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE))
  invisible(ok)
}

oc    <- readRDS(file.path(CFG$dir_out, "interim", "outcomes.rds"))
vis   <- readRDS(file.path(CFG$dir_out, "interim", "visits.rds"))
folds <- readRDS(file.path(CFG$dir_out, "interim", "folds.rds"))

d0 <- NULL; aud <- NULL; Pcv <- list()
section("0. Data and out-of-fold predictions", {
  d0  <- make_d0(oc, vis, folds, "in_set")
  aud <- add_groups(oc[oc$in_audit, ])
  cat("Modelling rows:", nrow(d0), "| audit rows:", nrow(aud), "\n")
  for (l in LMS) {
    dl <- d0[d0$landmark == l, ]
    Pcv[[l]] <- lapply(TIERS, function(f) run_cv(dl, f)$P)
  }
  cat("Out-of-fold predictions ready for", paste(LMS, collapse = ", "), "and tiers", paste(names(TIERS), collapse = ", "), "\n")
})

# ---- A1. bounds and tipping fraction ---------------------------------------
bounds <- NULL
section("A1. Observed risk bounds and tipping fraction (audit set, observable subjects)", {
  rows <- list()
  for (l in LMS) for (h in HS) {
    a  <- aud[aud$landmark == l, ]; cl <- a[[paste0("class_h", h)]]
    dl <- d0[d0$landmark == l, ];  dcl <- dl[[paste0("class_h", h)]]
    j  <- match(h, HORIZONS)
    for (v in c("all", GRP_VARS)) {
      lvA <- if (v == "all") rep("all", nrow(a)) else as.character(a[[v]])
      lvD <- if (v == "all") rep("all", nrow(dl)) else as.character(dl[[v]])
      for (u in setdiff(unique(lvA), c(NA, "Unknown"))) {
        ia <- which(lvA == u & cl != "administrative")
        if (length(ia) < MIN_N) next
        dem  <- sum(cl[ia] == "dementia"); lost <- sum(cl[ia] == "lost")
        lnv  <- sum(cl[ia] == "lost" & !a$has_followup[ia])
        idd  <- which(lvD == u)
        kmr  <- km_risk(dl$time[idd], dl$ev[idd], h)
        for (tn in names(TIERS)) {
          P  <- Pcv[[l]][[tn]][, j]
          io <- which(lvD == u & dcl != "administrative" & !is.na(P))
          pred <- mean(P[io])
          rows[[length(rows) + 1]] <- data.frame(landmark = l, horizon = h, variable = v, level = u, tier = tn,
            n_observable = length(ia), dementia = dem, lost = lost, lost_no_visit_after_landmark = lnv,
            risk_lower = dem / length(ia), risk_upper = (dem + lost) / length(ia),
            km_risk_modelling_set = kmr, mean_pred = pred,
            tipping_fraction = if (lost > 0) (pred * length(ia) - dem) / lost else NA_real_)
        }
      }
    }
  }
  bounds <- do.call(rbind, rows)
  bounds$read <- ifelse(is.na(bounds$tipping_fraction), "no lost subjects",
                 ifelse(bounds$tipping_fraction < 0, "observed already above predicted",
                 ifelse(bounds$tipping_fraction > 1, "no censoring scenario reconciles",
                        "censoring could reconcile")))
  out <- bounds; nn <- vapply(out, is.numeric, logical(1)); out[nn] <- lapply(out[nn], round, 3)
  save_csv(out, "06_bounds_tipping.csv")
  cat("Tier T0, 2-year horizon (tipping_fraction = share of lost subjects who must have converted for observed = predicted):\n")
  show(out[out$tier == "T0" & out$horizon == 2 & out$landmark == "bl",
           c("variable", "level", "n_observable", "dementia", "lost", "lost_no_visit_after_landmark",
             "risk_lower", "risk_upper", "mean_pred", "tipping_fraction", "read")], 40)
  cat("\nfull table: 06_bounds_tipping.csv (both tiers, bl and m12, 2 and 3 years)\n")
})

# ---- A2. covariate-adjusted IPCW -------------------------------------------
ipcw_tab <- NULL
section("A2. Covariate-adjusted IPCW observed risk vs Kaplan-Meier (modelling set)", {
  rows <- list()
  for (l in LMS) {
    dl <- d0[d0$landmark == l, ]
    dl$black <- as.integer(dl$race == "Black"); dl$hispanic <- as.integer(dl$hisp == "Hispanic")
    for (f in c("age", "male", "edu", "apoe4", "cdrsb")) dl[[f]][is.na(dl[[f]])] <- median(dl[[f]], na.rm = TRUE)
    fit_c <- survival::coxph(stats::as.formula(paste("survival::Surv(time, 1 - ev) ~", paste(CENS_COVARS, collapse = " + "))), data = dl)
    risk <- predict(fit_c, type = "risk")
    bh <- survival::basehaz(fit_c, centered = TRUE)
    H0 <- stepfun(bh$time, c(0, bh$hazard))
    GT <- pmax(exp(-H0(dl$time) * risk), 0.05)
    for (h in HS) {
      Gh <- pmax(exp(-H0(h) * risk), 0.05); j <- match(h, HORIZONS)
      for (v in c("all", GRP_VARS)) {
        lv <- if (v == "all") rep("all", nrow(dl)) else as.character(dl[[v]])
        for (u in setdiff(unique(lv), c(NA, "Unknown"))) {
          s <- which(lv == u)
          if (length(s) < 30) next
          kmr <- km_risk(dl$time[s], dl$ev[s], h)
          ipr <- ipcw_risk(dl$time, dl$ev, GT, Gh, h, idx = lv == u & !is.na(lv))
          for (tn in names(TIERS)) {
            pm <- mean(Pcv[[l]][[tn]][s, j], na.rm = TRUE)
            rows[[length(rows) + 1]] <- data.frame(landmark = l, horizon = h, variable = v, level = u, tier = tn,
              n = length(s), events_by_h = sum(dl$ev[s] == 1 & dl$time[s] <= h), km_risk = kmr, ipcw_risk = ipr,
              mean_pred = pm, OE_km = kmr / pm, OE_ipcw = ipr / pm)
          }
        }
      }
    }
  }
  ipcw_tab <- do.call(rbind, rows)
  out <- ipcw_tab; nn <- vapply(out, is.numeric, logical(1)); out[nn] <- lapply(out[nn], round, 3)
  save_csv(out, "06_ipcw_vs_km.csv")
  cat("Tier T0, baseline landmark, 2-year horizon:\n")
  show(out[out$tier == "T0" & out$landmark == "bl" & out$horizon == 2,
           c("variable", "level", "n", "events_by_h", "km_risk", "ipcw_risk", "mean_pred", "OE_km", "OE_ipcw")], 30)
  cat("\nIf OE_ipcw stays close to OE_km, covariate-explained censoring is not what drives the miscalibration.\n",
      "It cannot speak to censoring that depends on unmeasured factors; use the bounds above for that.\n", sep = "")
})

# ---- B. power --------------------------------------------------------------
section("B. Power: smallest O/E ratio detectable with 80% power (complete follow-up, 2-year horizon)", {
  power_cal <- function(P, m) {
    mu0 <- sum(P); s0 <- sqrt(sum(P * (1 - P)))
    pm <- pmin(1, m * P); mum <- sum(pm); sm <- sqrt(sum(pm * (1 - pm)))
    if (s0 == 0 || sm == 0) return(NA_real_)
    zc <- qnorm(0.975)
    pnorm((mu0 - zc * s0 - mum) / sm) + 1 - pnorm((mu0 + zc * s0 - mum) / sm)
  }
  rows <- list(); h <- 2; j <- match(h, HORIZONS)
  for (l in LMS) {
    dl <- d0[d0$landmark == l, ]; cl <- dl[[paste0("class_h", h)]]; P <- Pcv[[l]][["T0"]][, j]
    for (v in c("all", "sex", "race_bw", "hisp", "edu_grp", "apoe4_grp", "age_grp", "era")) {
      lv <- if (v == "all") rep("all", nrow(dl)) else as.character(dl[[v]])
      for (u in setdiff(unique(lv), c(NA, "Unknown", "missing"))) {
        s <- which(lv == u & cl %in% c("dementia", "death", "event_free") & !is.na(P))
        if (length(s) < MIN_N) next
        Ps <- P[s]
        lo <- seq(0.98, 0.05, by = -0.01); hi <- seq(1.02, 6, by = 0.02)
        plo <- sapply(lo, function(m) power_cal(Ps, m)); phi <- sapply(hi, function(m) power_cal(Ps, m))
        rows[[length(rows) + 1]] <- data.frame(landmark = l, variable = v, level = u,
          n_complete = length(s), expected_events = sum(Ps),
          min_detectable_OE_below = if (any(plo >= 0.8, na.rm = TRUE)) lo[which(plo >= 0.8)[1]] else NA_real_,
          min_detectable_OE_above = if (any(phi >= 0.8, na.rm = TRUE)) hi[which(phi >= 0.8)[1]] else NA_real_,
          power_at_0.7 = power_cal(Ps, 0.7), power_at_1.3 = power_cal(Ps, 1.3),
          type1_check = power_cal(Ps, 1))
      }
    }
  }
  pw <- do.call(rbind, rows)
  out <- pw; nn <- vapply(out, is.numeric, logical(1)); out[nn] <- lapply(out[nn], round, 2)
  save_csv(out, "06_power_by_subgroup.csv")
  show(out[out$landmark == "bl", setdiff(names(out), "type1_check")], 40)
  cat("\nmin_detectable_OE_below / above = closest O/E to 1 on each side that the group would flag as miscalibrated 80% of the time.\n",
      "NA means no O/E in the searched range (0.05 to 6) reaches 80% power. Ignores loss to follow-up, so it is optimistic.\n", sep = "")
  assign("pw", pw, envir = .GlobalEnv)
})

# ---- checks ----------------------------------------------------------------
section("C. Checks", {
  km_ok <- {
    dl <- d0[d0$landmark == "bl", ]; G <- cens_G(dl$time, dl$ev)
    abs(ipcw_risk(dl$time, dl$ev, G(dl$time), rep(G(2), nrow(dl)), 2) - km_risk(dl$time, dl$ev, 2)) < 0.01
  }
  check("IPCW with Kaplan-Meier censoring weights reproduces the Kaplan-Meier risk (within 0.01)", km_ok)
  check("risk bounds are ordered and inside 0 to 1", all(bounds$risk_lower <= bounds$risk_upper & bounds$risk_upper <= 1))
  check("the Kaplan-Meier risk lies between the bounds for the whole cohort at bl 2 years (within 0.05)", {
    b <- bounds[bounds$variable == "all" & bounds$landmark == "bl" & bounds$horizon == 2 & bounds$tier == "T0", ]
    b$km_risk_modelling_set >= b$risk_lower - 0.05 && b$km_risk_modelling_set <= b$risk_upper + 0.05 })
  check("power at O/E = 1 equals the 5% type-I rate", all(abs(pw$type1_check - 0.05) < 1e-6, na.rm = TRUE))
  check("IPCW risks are all between 0 and 1", all(ipcw_tab$ipcw_risk >= 0 & ipcw_tab$ipcw_risk <= 1, na.rm = TRUE))
  cat(sprintf("\n%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  save_csv(tests, "06_tests.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
