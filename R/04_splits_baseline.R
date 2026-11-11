# 04_splits_baseline.R -----------------------------------------------------
# 1. Subject-level 5-fold assignment shared by every model (stratified by
#    original protocol x whether the subject ever converts).
# 2. Leakage and pipeline tests, real and synthetic. Each prints PASS or FAIL.
# 3. Baseline landmark Cox model (age, APOE4 copies, CDR-SB), cross-validated.
#    Metrics per landmark and horizon: C-index, IPCW AUC, IPCW Brier, scaled
#    Brier, observed/expected ratio, plus O/E by subgroup.
#
# Death is treated as censoring here (time_cens / status_cens). Competing-risk
# versions come later. Missing predictors are imputed with TRAINING-fold medians.
#
# Writes: output/interim/folds.rds  (has subject IDs, stays local)
#         output/audit/04_*.csv, 04_baseline_report.txt   (counts and metrics only)
# Run:    source("R/04_splits_baseline.R")   (after 03_outcomes.R)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
if (!requireNamespace("survival", quietly = TRUE))
  stop("Package 'survival' is missing. It ships with R; run install.packages('survival').")

report <- file.path(CFG$dir_out, "audit", "04_baseline_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Splits and baseline | run:", format(Sys.time()), "\n")

# ---- settings -------------------------------------------------------------
K        <- 5
HORIZONS <- CFG$horizons_years
FEATS    <- c("age", "apoe4", "cdrsb")
MIN_N_SUBGROUP <- 30      # smallest subgroup reported in the O/E table
MIN_EV_SUBGROUP <- 5
# ---------------------------------------------------------------------------

tests <- data.frame(test = character(), result = character(), detail = character(),
                    stringsAsFactors = FALSE)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s %s\n", if (ok) "PASS" else "FAIL", name, detail))
  tests <<- rbind(tests, data.frame(test = name, result = if (ok) "PASS" else "FAIL",
                                    detail = detail, stringsAsFactors = FALSE))
  invisible(ok)
}

oc  <- readRDS(file.path(CFG$dir_out, "interim", "outcomes.rds"))
vis <- readRDS(file.path(CFG$dir_out, "interim", "visits.rds"))

# ---- 1. analysis table ----------------------------------------------------
d0 <- NULL
section("1. Analysis table (modelling set, one row per subject and landmark)", {
  d0 <- oc %>% filter(in_set) %>%
    left_join(vis %>% select(RID, EXAMDATE, cdrsb, apoe4), by = c("RID", "t0" = "EXAMDATE")) %>%
    mutate(time = time_cens, ev = status_cens, age = age_lm,
           race_bw = ifelse(race %in% c("Black", "White"), race, NA_character_))
  cat("Rows:", nrow(d0), "| subjects:", n_distinct(d0$RID), "\n")
  miss <- d0 %>% group_by(landmark) %>%
    summarise(n = n(), events_any = sum(ev), missing_cdrsb = sum(is.na(cdrsb)),
              missing_apoe4 = sum(is.na(apoe4)), missing_age = sum(is.na(age))) %>% ungroup()
  show(miss); save_csv(miss, "04_predictor_missingness.csv")
})

# ---- 2. folds --------------------------------------------------------------
folds <- NULL
section("2. Subject-level folds", {
  subj <- d0 %>% group_by(RID) %>%
    summarise(ever_dem = any(ev == 1), origprot = first(origprot)) %>% ungroup() %>% arrange(RID)
  subj$stratum <- paste(subj$origprot, ifelse(subj$ever_dem, "converter", "non-converter"))
  set.seed(CFG$seed)
  subj$fold <- NA_integer_
  counter <- 0L
  for (s in sort(unique(subj$stratum))) {
    idx <- which(subj$stratum == s)
    idx <- idx[sample.int(length(idx))]
    subj$fold[idx] <- ((counter + seq_along(idx) - 1L) %% K) + 1L
    counter <- counter + length(idx)
  }
  folds <- subj[, c("RID", "fold", "stratum")]
  saveRDS(folds, file.path(CFG$dir_out, "interim", "folds.rds"))
  d0 <- d0 %>% left_join(folds[, c("RID", "fold")], by = "RID")
  bal <- d0 %>% group_by(fold) %>%
    summarise(n_subjects = n_distinct(RID), n_rows = n(),
              n_bl = sum(landmark == "bl"), converters_bl = sum(landmark == "bl" & ev == 1),
              n_black_subjects = n_distinct(RID[race %in% "Black"]),
              n_hispanic_subjects = n_distinct(RID[hisp %in% "Hispanic"])) %>% ungroup()
  show(bal); save_csv(bal, "04_fold_balance.csv")
  cat("\nSeed:", CFG$seed, "| K =", K, "| strata:", length(unique(subj$stratum)), "\n")
})

# ---- 3. structural leakage tests (real data) -------------------------------
section("3. Leakage and pipeline tests on the real tables", {
  cat("Fold assignment\n")
  disjoint <- function(df) {
    sets <- split(df$RID, df$fold)
    ks <- names(sets)
    for (a in seq_along(ks)) for (b in seq_along(ks)) if (a < b)
      if (length(intersect(sets[[a]], sets[[b]])) > 0) return(FALSE)
    TRUE
  }
  check("every subject sits in exactly one fold", all(tapply(d0$fold, d0$RID, function(x) length(unique(x))) == 1))
  check("fold subject sets are pairwise disjoint (analysis table)", disjoint(unique(d0[, c("RID", "fold")])))
  check("no row without a fold", !any(is.na(d0$fold)))
  sz <- table(folds$fold)
  check("fold sizes within 10% of each other", max(sz) / min(sz) <= 1.10,
        paste0("(", paste(sz, collapse = "/"), ")"))
  bad <- data.frame(RID = rep(1:50, each = 4), fold = sample(1:K, 200, replace = TRUE))
  check("the disjointness check catches a record-level split (synthetic)", !disjoint(bad))

  cat("\nTime ordering\n")
  check("no landmark row has a dementia diagnosis on or before time zero",
        all(is.na(d0$t_dem) | d0$t_dem > d0$t0))
  check("all follow-up times are positive", all(d0$time > 0))
  vd <- vis %>% select(RID, EXAMDATE, DIAGNOSIS)
  chk <- d0 %>% select(RID, t0) %>% left_join(vd, by = c("RID", "t0" = "EXAMDATE"))
  check("every landmark row is an MCI diagnosis at its landmark date", all(chk$DIAGNOSIS == 2, na.rm = TRUE) && !anyNA(chk$DIAGNOSIS))
  dcols <- paste0(c("cdr", "adas", "faq", "mmse", "ravlt", "fs", "plasma"), "_date")
  dcols <- dcols[dcols %in% names(vis)]
  late <- sapply(dcols, function(cn) sum(vis[[cn]] > vis$EXAMDATE, na.rm = TRUE))
  check("no assessment in the visit table is dated after its visit", all(late == 0),
        paste0("(checked ", paste(dcols, collapse = ", "), ")"))
  lm_rows <- d0 %>% select(RID, t0) %>% left_join(vis, by = c("RID", "t0" = "EXAMDATE"))
  late_lm <- sapply(dcols, function(cn) sum(lm_rows[[cn]] > lm_rows$t0, na.rm = TRUE))
  check("no feature at any landmark row is dated after time zero", all(late_lm == 0))
})

# ---- baseline machinery ----------------------------------------------------
run_cv <- function(d, feats) {
  n <- nrow(d); P <- matrix(NA_real_, n, length(HORIZONS)); fc <- rep(NA_real_, K)
  fml <- stats::as.formula(paste("survival::Surv(time, ev) ~", paste(feats, collapse = " + ")))
  for (k in seq_len(K)) {
    tr <- d$fold != k; te <- !tr
    if (sum(d$ev[tr]) < 5 || sum(te) < 2) next
    med <- sapply(feats, function(f) median(d[[f]][tr], na.rm = TRUE))
    imp <- function(x) { for (f in feats) x[[f]][is.na(x[[f]])] <- med[[f]]; x }
    dtr <- imp(d[tr, ]); dte <- imp(d[te, ])
    fit <- survival::coxph(fml, data = dtr)
    lp  <- predict(fit, newdata = dte, type = "lp")
    sf  <- survival::survfit(fit, newdata = dte)
    for (j in seq_along(HORIZONS))
      P[te, j] <- 1 - as.numeric(summary(sf, times = HORIZONS[j], extend = TRUE)$surv)
    fc[k] <- survival::concordance(survival::Surv(time, ev) ~ lp,
                                   data = data.frame(time = dte$time, ev = dte$ev, lp = lp),
                                   reverse = TRUE)$concordance
  }
  list(P = P, fold_c = fc)
}

cens_G <- function(time, ev) {
  km <- survival::survfit(survival::Surv(time, 1 - ev) ~ 1)
  if (length(km$time) == 0) return(function(t) rep(1, length(t)))
  sf <- stepfun(km$time, c(1, km$surv))
  function(t) pmax(sf(t), 0.05)
}
km_risk <- function(time, ev, h) {
  if (sum(ev) == 0) return(0)
  1 - summary(survival::survfit(survival::Surv(time, ev) ~ 1), times = h, extend = TRUE)$surv
}
metrics_h <- function(time, ev, p, h) {
  ok <- !is.na(p); time <- time[ok]; ev <- ev[ok]; p <- p[ok]; n <- length(p)
  G <- cens_G(time, ev)
  case <- ev == 1 & time <= h; ctrl <- time > h
  w <- numeric(n); w[case] <- 1 / G(time[case]); w[ctrl] <- 1 / G(h)
  y <- as.numeric(case)
  bs <- sum(w * (y - p)^2) / n
  p0 <- km_risk(time, ev, h)
  bs0 <- sum(w * (y - p0)^2) / n
  auc <- NA_real_
  if (any(case) && any(ctrl)) {
    pc <- p[case]; pk <- p[ctrl]; wc <- w[case]; wk <- w[ctrl]
    num <- sum(outer(wc, wk) * ((outer(pc, pk, ">")) + 0.5 * outer(pc, pk, "==")))
    auc <- num / (sum(wc) * sum(wk))
  }
  data.frame(n = n, events_by_h = sum(case), AUC = auc, Brier = bs, Brier_null = bs0,
             scaled_Brier = 1 - bs / bs0, mean_pred = mean(p), km_observed = p0,
             O_over_E = p0 / mean(p))
}

# ---- 4. baseline Cox, cross-validated -------------------------------------
oof <- list()
section("4. Baseline landmark Cox (age, APOE4 copies, CDR-SB), 5-fold subject-level CV", {
  res <- list()
  for (l in CFG$landmarks) {
    dl <- d0[d0$landmark == l, ]
    cv <- run_cv(dl, FEATS)
    oof[[l]] <- list(d = dl, P = cv$P)
    for (j in seq_along(HORIZONS)) {
      m <- metrics_h(dl$time, dl$ev, cv$P[, j], HORIZONS[j])
      m$C_mean_folds <- mean(cv$fold_c, na.rm = TRUE); m$C_sd_folds <- sd(cv$fold_c, na.rm = TRUE)
      res[[length(res) + 1]] <- cbind(landmark = l, horizon = HORIZONS[j], m)
    }
  }
  out <- do.call(rbind, res)
  num <- vapply(out, is.numeric, logical(1)); out[num] <- lapply(out[num], round, 3)
  show(out, 40); save_csv(out, "04_baseline_cox_metrics.csv")
  cat("\nHeadline numbers are out-of-fold. C-index is the mean over folds. AUC and Brier use IPCW with the censoring\n",
      "distribution estimated on the pooled out-of-fold rows. scaled_Brier = 1 - Brier / Brier of a constant (KM) risk.\n", sep = "")

  hr <- do.call(rbind, lapply(CFG$landmarks, function(l) {
    dl <- d0[d0$landmark == l, ]
    for (f in FEATS) dl[[f]][is.na(dl[[f]])] <- median(dl[[f]], na.rm = TRUE)
    fit <- survival::coxph(stats::as.formula(paste("survival::Surv(time, ev) ~", paste(FEATS, collapse = " + "))), data = dl)
    ci <- summary(fit)$conf.int
    data.frame(landmark = l, term = rownames(ci), HR = round(ci[, 1], 3),
               lcl = round(ci[, 3], 3), ucl = round(ci[, 4], 3), row.names = NULL)
  }))
  cat("\nHazard ratios (full-data fit, for interpretation only; age per year, CDR-SB per point, APOE4 per copy):\n")
  show(hr, 20); save_csv(hr, "04_baseline_hazard_ratios.csv")
})

# ---- 5. negative control ---------------------------------------------------
section("5. Negative control: outcomes shuffled across subjects", {
  set.seed(CFG$seed + 1)
  for (l in CFG$landmarks) {
    dl <- d0[d0$landmark == l, ]
    idx <- sample.int(nrow(dl)); dl$time <- dl$time[idx]; dl$ev <- dl$ev[idx]
    cv <- run_cv(dl, FEATS)
    cm <- mean(cv$fold_c, na.rm = TRUE)
    tol <- if (l %in% c("bl", "m12")) 0.06 else 0.10
    check(paste0("shuffled outcomes give C near 0.5 at ", l), abs(cm - 0.5) < tol,
          sprintf("(C = %.3f, tolerance +/- %.2f)", cm, tol))
  }
})

# ---- 6. observed / expected by subgroup ------------------------------------
section("6. Calibration by subgroup (O/E = KM-observed risk / mean predicted risk)", {
  vars <- c("sex", "race_bw", "hisp", "edu_grp", "apoe4_grp", "age_grp", "origprot")
  rows <- list()
  for (l in CFG$landmarks) {
    dl <- oof[[l]]$d; P <- oof[[l]]$P
    for (j in seq_along(HORIZONS)) for (v in vars) {
      lv <- as.character(dl[[v]])
      for (u in setdiff(unique(lv), c(NA, "Unknown", "missing"))) {
        s <- which(lv == u & !is.na(P[, j]))
        if (length(s) < MIN_N_SUBGROUP) next
        h <- HORIZONS[j]
        obs <- km_risk(dl$time[s], dl$ev[s], h); ex <- mean(P[s, j])
        ne <- sum(dl$ev[s] == 1 & dl$time[s] <= h)
        if (ne < MIN_EV_SUBGROUP) next
        rows[[length(rows) + 1]] <- data.frame(landmark = l, horizon = h, variable = v, level = u,
          n = length(s), events_by_h = ne, km_observed = round(obs, 3), mean_pred = round(ex, 3),
          O_over_E = round(obs / ex, 2))
      }
    }
  }
  oe <- do.call(rbind, rows)
  show(oe %>% filter(landmark == "bl", horizon == 2), 40)
  save_csv(oe, "04_baseline_oe_by_subgroup.csv")
  cat("\nSubgroups under", MIN_N_SUBGROUP, "subjects or", MIN_EV_SUBGROUP, "events are left out. Intervals are not shown yet; read O/E for small groups with care.\n")
})

# ---- 7. synthetic positive control -----------------------------------------
section("7. Synthetic positive control", {
  set.seed(CFG$seed)
  n <- 3000; x <- rnorm(n)
  tt <- rexp(n, rate = 0.1 * exp(0.7 * x)); cc <- runif(n, 0, 15)
  time <- pmin(tt, cc); ev <- as.integer(tt <= cc)
  fit <- survival::coxph(survival::Surv(time, ev) ~ x)
  b <- unname(coef(fit)); cidx <- survival::concordance(fit)$concordance
  check("Cox recovers a known log hazard ratio of 0.7 (+/- 0.10)", abs(b - 0.7) < 0.10,
        sprintf("(estimate %.3f)", b))
  check("C-index on the synthetic data is above 0.65", cidx > 0.65, sprintf("(C = %.3f)", cidx))
  pr <- 1 - exp(-0.1 * exp(0.7 * x) * 5)
  mm <- metrics_h(time, ev, pr, 5)
  check("IPCW O/E is near 1 when predictions are the true risks (0.9 to 1.1)",
        mm$O_over_E > 0.9 && mm$O_over_E < 1.1, sprintf("(O/E = %.2f)", mm$O_over_E))
})

# ---- 8. test summary -------------------------------------------------------
section("8. Test summary", {
  cat(sprintf("%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  if (any(tests$result == "FAIL")) { cat("FAILED:\n"); print(tests[tests$result == "FAIL", ], row.names = FALSE) }
  save_csv(tests, "04_tests.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
