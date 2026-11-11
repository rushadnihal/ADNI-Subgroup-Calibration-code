# 05_tiers_bootstrap.R -----------------------------------------------------
# 1. Feature tiers, cross-validated with the SAME folds as 04:
#      T0   age, APOE4, CDR-SB                      (the 04 baseline, as a regression check)
#      A    age, sex, education, APOE4
#      B    A + CDR-SB, ADAS-Cog13, RAVLT immediate recall, FAQ
#      C    B + hippocampus, entorhinal, middle temporal, ventricles (each per 1000 of ICV)
#      T0+era, B+era   the same with a protocol-era term (ADNI1 / ADNIGO+ADNI2 / ADNI3+ADNI4),
#                      bl and m12 only, to test whether protocol miscalibration is a case-mix effect
# 2. Subject-level bootstrap intervals (B resamples of the out-of-fold predictions) for AUC,
#    scaled Brier and O/E, for AUC gains over T0, and for subgroup O/E.
#
# The bootstrap resamples evaluation rows (one row per subject within a landmark). It does not
# refit the models, so it captures evaluation noise, not training noise.
# Death is treated as censoring (time_cens / status_cens), as in 04.
#
# Writes: output/audit/05_*.csv, 05_tiers_report.txt (counts and metrics only)
# Run:    source("R/05_tiers_bootstrap.R")   (after 04). About 4 to 8 minutes; set B lower to go faster.
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
if (!requireNamespace("survival", quietly = TRUE))
  stop("Package 'survival' is missing. It ships with R; run install.packages('survival').")

report <- file.path(CFG$dir_out, "audit", "05_tiers_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Tiers and bootstrap | run:", format(Sys.time()), "\n")

# ---- settings -------------------------------------------------------------
B        <- 300
K        <- 5
HORIZONS <- CFG$horizons_years
TIERS <- list(
  "T0" = c("age", "apoe4", "cdrsb"),
  "A"  = c("age", "male", "edu", "apoe4"),
  "B"  = c("age", "male", "edu", "apoe4", "cdrsb", "adas13", "ravlt_immediate", "faq"),
  "C"  = c("age", "male", "edu", "apoe4", "cdrsb", "adas13", "ravlt_immediate", "faq",
           "hippo", "entorhinal", "midtemp", "ventricles"),
  "T0+era" = c("age", "apoe4", "cdrsb", "era"),
  "B+era"  = c("age", "male", "edu", "apoe4", "cdrsb", "adas13", "ravlt_immediate", "faq", "era"))
ERA_TIERS <- c("T0+era", "B+era"); ERA_LM <- c("bl", "m12")
SG_TIERS  <- c("T0", "B", "T0+era", "B+era"); SG_LM <- c("bl", "m12"); SG_H <- 2
SG_VARS   <- c("sex", "race_bw", "hisp", "edu_grp", "apoe4_grp", "age_grp", "origprot")
MIN_N_SUBGROUP <- 30; MIN_EV_SUBGROUP <- 5
# ---------------------------------------------------------------------------

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

# ---- evaluation machinery (same definitions as 04) -------------------------
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
met3 <- function(time, ev, p, h) {
  ok <- !is.na(p)
  if (!any(ok)) return(c(AUC = NA_real_, sBrier = NA_real_, OE = NA_real_))
  time <- time[ok]; ev <- ev[ok]; p <- p[ok]; n <- length(p)
  G <- cens_G(time, ev)
  case <- ev == 1 & time <= h; ctrl <- time > h
  w <- numeric(n); w[case] <- 1 / G(time[case]); w[ctrl] <- 1 / G(h)
  y <- as.numeric(case)
  p0 <- km_risk(time, ev, h)
  bs <- sum(w * (y - p)^2) / n; bs0 <- sum(w * (y - p0)^2) / n
  auc <- NA_real_
  if (any(case) && any(ctrl)) {
    pc <- p[case]; pk <- p[ctrl]; wc <- w[case]; wk <- w[ctrl]
    auc <- sum(outer(wc, wk) * ((outer(pc, pk, ">")) + 0.5 * outer(pc, pk, "=="))) / (sum(wc) * sum(wk))
  }
  c(AUC = auc, sBrier = 1 - bs / bs0, OE = p0 / mean(p))
}
run_cv <- function(d, feats) {
  n <- nrow(d); P <- matrix(NA_real_, n, length(HORIZONS)); fc <- rep(NA_real_, K)
  fml <- stats::as.formula(paste("survival::Surv(time, ev) ~", paste(feats, collapse = " + ")))
  for (k in seq_len(K)) {
    tr <- d$fold != k; te <- !tr
    if (sum(d$ev[tr]) < 5 || sum(te) < 2) next
    num <- feats[sapply(feats, function(f) is.numeric(d[[f]]))]
    med <- sapply(num, function(f) median(d[[f]][tr], na.rm = TRUE))
    imp <- function(x) { for (f in num) x[[f]][is.na(x[[f]])] <- med[[f]]; x }
    dtr <- imp(d[tr, ]); dte <- imp(d[te, ])
    ok <- tryCatch({
      fit <- survival::coxph(fml, data = dtr)
      lp  <- predict(fit, newdata = dte, type = "lp")
      sf  <- survival::survfit(fit, newdata = dte)
      for (j in seq_along(HORIZONS))
        P[te, j] <- 1 - as.numeric(summary(sf, times = HORIZONS[j], extend = TRUE)$surv)
      fc[k] <- survival::concordance(survival::Surv(time, ev) ~ lp,
                 data = data.frame(time = dte$time, ev = dte$ev, lp = lp), reverse = TRUE)$concordance
      TRUE
    }, error = function(e) FALSE)
    if (!ok) cat("  note: model failed in fold", k, "\n")
  }
  list(P = P, fold_c = fc)
}
q95 <- function(x) { x <- x[is.finite(x)]; if (length(x) < 20) c(NA_real_, NA_real_) else unname(quantile(x, c(0.025, 0.975))) }

# ---- 1. analysis table -----------------------------------------------------
d0 <- NULL
section("1. Analysis table with tier predictors", {
  vv <- vis %>% select(RID, EXAMDATE, cdrsb, apoe4, edu, adas13, ravlt_immediate, faq,
                       hippocampus, entorhinal, mid_temporal, ventricles, icv)
  d0 <- oc %>% filter(in_set) %>%
    left_join(vv, by = c("RID", "t0" = "EXAMDATE")) %>%
    left_join(folds[, c("RID", "fold")], by = "RID") %>%
    mutate(time = time_cens, ev = status_cens, age = age_lm,
           male = as.integer(sex == "Male"),
           hippo = hippocampus / icv * 1000, midtemp = mid_temporal / icv * 1000,
           entorhinal = entorhinal / icv * 1000, ventricles = ventricles / icv * 1000,
           race_bw = ifelse(race %in% c("Black", "White"), race, NA_character_),
           era = factor(case_when(origprot == "ADNI1" ~ "ADNI1",
                                  origprot %in% c("ADNIGO", "ADNI2") ~ "ADNIGO/2",
                                  TRUE ~ "ADNI3/4"), levels = c("ADNI1", "ADNIGO/2", "ADNI3/4")))
  cat("Rows:", nrow(d0), "| subjects:", n_distinct(d0$RID), "| rows without a fold:", sum(is.na(d0$fold)), "\n")
  allf <- unique(unlist(TIERS)); allf <- setdiff(allf, "era")
  ms <- do.call(rbind, lapply(CFG$landmarks, function(l) {
    x <- d0[d0$landmark == l, ]
    data.frame(landmark = l, n = nrow(x), events_any = sum(x$ev),
               t(sapply(allf, function(f) round(100 * mean(is.na(x[[f]])), 1))), check.names = FALSE)
  }))
  cat("\nPercent missing at each landmark (imputed with training-fold medians):\n")
  show(ms, 10); save_csv(ms, "05_feature_missingness.csv")
  epv <- data.frame(tier = names(TIERS), n_terms = sapply(TIERS, length))
  epv$events_per_term_m36 <- round(sum(d0$ev[d0$landmark == "m36"]) / epv$n_terms, 1)
  cat("\nEvents per model term at m36 (the thinnest landmark; below ~10 means overfitting risk):\n"); show(epv, 10)
})

# ---- 2. cross-validate every tier, bootstrap, assemble ----------------------
main <- list(); delta <- list(); sg <- list(); fcl <- list()
section("2. Tiers, out-of-fold metrics and bootstrap intervals", {
  for (l in CFG$landmarks) {
    dl <- d0[d0$landmark == l, ]; n <- nrow(dl)
    Ps <- list()
    for (tn in names(TIERS)) {
      if (tn %in% ERA_TIERS && !(l %in% ERA_LM)) next
      cv <- run_cv(dl, TIERS[[tn]]); Ps[[tn]] <- cv$P
      fcl[[paste(l, tn)]] <- data.frame(landmark = l, tier = tn, C_mean_folds = mean(cv$fold_c, na.rm = TRUE),
                                        C_sd_folds = sd(cv$fold_c, na.rm = TRUE))
    }
    tn_l <- names(Ps)
    bt <- array(NA_real_, c(B, length(tn_l), length(HORIZONS), 3),
                dimnames = list(NULL, tn_l, as.character(HORIZONS), c("AUC", "sBrier", "OE")))
    # subgroup levels (defined once on the original sample)
    sg_levels <- list()
    if (l %in% SG_LM) {
      j2 <- which(HORIZONS == SG_H)
      for (v in SG_VARS) {
        lv <- as.character(dl[[v]])
        for (u in setdiff(unique(lv), c(NA, "Unknown", "missing"))) {
          s <- which(lv == u)
          if (length(s) >= MIN_N_SUBGROUP && sum(dl$ev[s] == 1 & dl$time[s] <= SG_H) >= MIN_EV_SUBGROUP)
            sg_levels[[length(sg_levels) + 1]] <- list(var = v, level = u, lv = lv)
        }
      }
    }
    sg_bt <- lapply(intersect(SG_TIERS, tn_l), function(tn) matrix(NA_real_, B, length(sg_levels)))
    names(sg_bt) <- intersect(SG_TIERS, tn_l)
    set.seed(CFG$seed + 100 + match(l, CFG$landmarks))
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      for (tn in tn_l) for (j in seq_along(HORIZONS))
        bt[b, tn, j, ] <- met3(dl$time[idx], dl$ev[idx], Ps[[tn]][idx, j], HORIZONS[j])
      if (length(sg_levels)) for (tn in names(sg_bt)) for (m in seq_along(sg_levels)) {
        u <- sg_levels[[m]]; s <- idx[which(u$lv[idx] == u$level)]
        if (length(s) >= 5) sg_bt[[tn]][b, m] <- km_risk(dl$time[s], dl$ev[s], SG_H) / mean(Ps[[tn]][s, j2], na.rm = TRUE)
      }
      if (b %% 100 == 0) cat("  ", l, ": bootstrap", b, "of", B, "\n")
    }
    for (tn in tn_l) for (j in seq_along(HORIZONS)) {
      pt <- met3(dl$time, dl$ev, Ps[[tn]][, j], HORIZONS[j])
      ci <- lapply(c("AUC", "sBrier", "OE"), function(k) q95(bt[, tn, j, k]))
      main[[length(main) + 1]] <- data.frame(landmark = l, tier = tn, horizon = HORIZONS[j], n = n,
        events_by_h = sum(dl$ev == 1 & dl$time <= HORIZONS[j]), n_terms = length(TIERS[[tn]]),
        AUC = pt[["AUC"]], AUC_lcl = ci[[1]][1], AUC_ucl = ci[[1]][2],
        scaled_Brier = pt[["sBrier"]], sBrier_lcl = ci[[2]][1], sBrier_ucl = ci[[2]][2],
        O_over_E = pt[["OE"]], OE_lcl = ci[[3]][1], OE_ucl = ci[[3]][2])
      if (tn != "T0") {
        dd <- bt[, tn, j, "AUC"] - bt[, "T0", j, "AUC"]
        p0 <- met3(dl$time, dl$ev, Ps[["T0"]][, j], HORIZONS[j])[["AUC"]]
        cd <- q95(dd)
        delta[[length(delta) + 1]] <- data.frame(landmark = l, tier = tn, horizon = HORIZONS[j],
          delta_AUC_vs_T0 = pt[["AUC"]] - p0, lcl = cd[1], ucl = cd[2])
      }
    }
    for (tn in names(sg_bt)) for (m in seq_along(sg_levels)) {
      u <- sg_levels[[m]]; s <- which(u$lv == u$level)
      ci <- q95(sg_bt[[tn]][, m])
      sg[[length(sg) + 1]] <- data.frame(landmark = l, tier = tn, horizon = SG_H, variable = u$var, level = u$level,
        n = length(s), events_by_h = sum(dl$ev[s] == 1 & dl$time[s] <= SG_H),
        km_observed = km_risk(dl$time[s], dl$ev[s], SG_H), mean_pred = mean(Ps[[tn]][s, j2], na.rm = TRUE),
        O_over_E = km_risk(dl$time[s], dl$ev[s], SG_H) / mean(Ps[[tn]][s, j2], na.rm = TRUE),
        lcl = ci[1], ucl = ci[2])
    }
  }
  rd <- function(x) { nn <- vapply(x, is.numeric, logical(1)); x[nn] <- lapply(x[nn], round, 3); x }
  main_df  <- rd(do.call(rbind, main));  delta_df <- rd(do.call(rbind, delta))
  sg_df    <- rd(do.call(rbind, sg));    fc_df    <- rd(do.call(rbind, fcl))
  main_df  <- merge(main_df, fc_df, by = c("landmark", "tier"), sort = FALSE)
  main_df  <- main_df[order(match(main_df$landmark, CFG$landmarks), match(main_df$tier, names(TIERS)), main_df$horizon), ]
  save_csv(main_df, "05_tier_metrics.csv"); save_csv(delta_df, "05_auc_gain_vs_T0.csv"); save_csv(sg_df, "05_subgroup_oe_ci.csv")
  assign("main_df", main_df, envir = .GlobalEnv)
  cat("\nTier metrics, 2-year horizon (95% bootstrap intervals; full table in 05_tier_metrics.csv):\n")
  show(main_df[main_df$horizon == 2, c("landmark", "tier", "n_terms", "events_by_h", "AUC", "AUC_lcl", "AUC_ucl",
                                       "scaled_Brier", "O_over_E", "C_mean_folds")], 40)
  cat("\nAUC gain over T0 (paired bootstrap), 2-year horizon:\n")
  show(delta_df[delta_df$horizon == 2, ], 30)
  cat("\nSubgroup O/E at baseline, 2-year horizon, with intervals:\n")
  show(sg_df[sg_df$landmark == "bl", c("tier", "variable", "level", "n", "events_by_h", "km_observed", "mean_pred", "O_over_E", "lcl", "ucl")], 80)
})

# ---- 3. checks -------------------------------------------------------------
section("3. Checks", {
  ref <- tryCatch(readr::read_csv(file.path(CFG$dir_out, "audit", "04_baseline_cox_metrics.csv"),
                                  show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(ref)) {
    m <- merge(main_df[main_df$tier == "T0", c("landmark", "horizon", "AUC")],
               as.data.frame(ref[, c("landmark", "horizon", "AUC")]), by = c("landmark", "horizon"), suffixes = c("_05", "_04"))
    check("T0 AUC reproduces the 04 baseline at every landmark and horizon (within 0.002)",
          nrow(m) == 12 && all(abs(m$AUC_05 - m$AUC_04) < 0.002),
          sprintf("(largest difference %.4f)", max(abs(m$AUC_05 - m$AUC_04))))
  } else check("04 baseline file found for the regression check", FALSE)
  check("no tier has an AUC below 0.5 at any landmark and horizon", all(main_df$AUC > 0.5, na.rm = TRUE))
  check("bootstrap intervals contain their point estimates (AUC)",
        all(main_df$AUC >= main_df$AUC_lcl - 0.02 & main_df$AUC <= main_df$AUC_ucl + 0.02, na.rm = TRUE))
  check("no row without a fold", !any(is.na(d0$fold)))
  cat(sprintf("\n%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  save_csv(tests, "05_tests.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
