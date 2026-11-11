# 07_model_classes.R -------------------------------------------------------
# Benchmark of model classes on one cohort, one set of subject-level folds, one
# feature set and one metric set (landmarks bl, m12, m24, m36; horizons 1, 2, 3 y).
#
#   Cox-B          Cox, Tier B (age, sex, education, APOE4, CDR-SB, ADAS-Cog13, RAVLT, FAQ)   [reference]
#   Cox-C          Cox, Tier C (B + hippocampus, entorhinal, middle temporal, ventricles / ICV)
#   GLMNET-C       elastic-net Cox (alpha 0.5, lambda by inner 5-fold CV on the training fold), Tier C
#   RSF-C          random survival forest (ranger, 500 trees, min node size 15, untuned), Tier C
#   ... +era       the same with a protocol-era term (bl and m12 only)
#
# Every model is trained only on the training folds; median imputation uses training-fold medians.
# Metrics: IPCW AUC, scaled Brier, O/E, each with subject-level bootstrap intervals over the
# out-of-fold predictions (models are not refit), AUC gain over Cox-B, and subgroup O/E at
# bl and m12 (2-year horizon). Death is treated as censoring. Hyperparameters are not tuned.
#
# Writes: output/interim/oof_predictions.rds (has subject IDs, stays local)
#         output/audit/07_*.csv, 07_model_classes_report.txt
# Run:    source("R/07_model_classes.R")   (after 06). About 8 to 15 minutes; lower B to go faster.
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
for (p in c("survival", "glmnet", "ranger"))
  if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is missing; run install.packages('", p, "').")
K <- 5; HORIZONS <- CFG$horizons_years
source(if (file.exists("R/utils_eval.R")) "R/utils_eval.R" else "utils_eval.R")

report <- file.path(CFG$dir_out, "audit", "07_model_classes_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Model classes | run:", format(Sys.time()), "\n")

# ---- settings -------------------------------------------------------------
B  <- 300
FB <- c("age", "male", "edu", "apoe4", "cdrsb", "adas13", "ravlt_immediate", "faq")
FC <- c(FB, "hippo", "entorhinal", "midtemp", "ventricles")
MODELS <- list(
  list(name = "Cox-B",        feats = FB,              type = "cox"),
  list(name = "Cox-C",        feats = FC,              type = "cox"),
  list(name = "GLMNET-C",     feats = FC,              type = "glmnet"),
  list(name = "RSF-C",        feats = FC,              type = "rsf"),
  list(name = "Cox-B+era",    feats = c(FB, "era"),    type = "cox",    era = TRUE),
  list(name = "GLMNET-C+era", feats = c(FC, "era"),    type = "glmnet", era = TRUE),
  list(name = "RSF-C+era",    feats = c(FC, "era"),    type = "rsf",    era = TRUE))
REF <- "Cox-B"; ERA_LM <- c("bl", "m12")
SG_LM <- c("bl", "m12"); SG_H <- 2
SG_VARS <- c("sex", "race_bw", "hisp", "edu_grp", "apoe4_grp", "age_grp", "origprot")
MIN_N_SUBGROUP <- 30; MIN_EV_SUBGROUP <- 5
# ---------------------------------------------------------------------------

tests <- data.frame(test = character(), result = character(), detail = character(), stringsAsFactors = FALSE)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s %s\n", if (ok) "PASS" else "FAIL", name, detail))
  tests <<- rbind(tests, data.frame(test = name, result = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE))
  invisible(ok)
}

# ---- learners ----------------------------------------------------------------
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
learn_rsf <- function(dtr, dte, feats) {
  tr <- dtr[, c("time", "ev", feats)]; te <- dte[, feats, drop = FALSE]
  if ("era" %in% feats) { tr$era <- as.integer(tr$era); te$era <- as.integer(te$era) }
  fit <- ranger::ranger(dependent.variable.name = "time", status.variable.name = "ev", data = tr,
                        num.trees = 500, min.node.size = 15, seed = CFG$seed, num.threads = 1)
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
predict_model <- function(d, m) {
  switch(m$type,
         cox    = run_cv(d, m$feats)$P,
         glmnet = run_cv_gen(d, m$feats, learn_glmnet),
         rsf    = run_cv_gen(d, m$feats, learn_rsf))
}

oc    <- readRDS(file.path(CFG$dir_out, "interim", "outcomes.rds"))
vis   <- readRDS(file.path(CFG$dir_out, "interim", "visits.rds"))
folds <- readRDS(file.path(CFG$dir_out, "interim", "folds.rds"))

d0 <- NULL; Pall <- list()
section("1. Out-of-fold predictions for every model class", {
  d0 <- make_d0(oc, vis, folds, "in_set")
  cat("Rows:", nrow(d0), "| subjects:", length(unique(d0$RID)), "| rows without a fold:", sum(is.na(d0$fold)), "\n")
  for (l in CFG$landmarks) {
    dl <- d0[d0$landmark == l, ]; Pall[[l]] <- list()
    for (m in MODELS) {
      if (isTRUE(m$era) && !(l %in% ERA_LM)) next
      t0 <- Sys.time()
      Pall[[l]][[m$name]] <- predict_model(dl, m)
      cat(sprintf("  %-4s %-13s %5.1f s\n", l, m$name, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }
})

# ---- glmnet coefficients (full-data fit, for interpretation) -----------------
section("2. Elastic-net Cox: terms kept at each landmark (full-data fit, per 1 SD)", {
  rows <- list()
  for (l in CFG$landmarks) {
    dl <- d0[d0$landmark == l, ]
    for (f in FC) dl[[f]][is.na(dl[[f]])] <- median(dl[[f]], na.rm = TRUE)
    x <- as.matrix(dl[, FC]); y <- cbind(time = dl$time, status = dl$ev)
    set.seed(CFG$seed)
    cvf <- glmnet::cv.glmnet(x, y, family = "cox", alpha = 0.5, nfolds = 5)
    cf <- as.numeric(coef(cvf, s = "lambda.min")); sdx <- apply(x, 2, sd)
    rows[[l]] <- data.frame(landmark = l, term = FC, HR_per_SD = round(exp(cf * sdx), 3), kept = cf != 0,
                            lambda_min = signif(cvf$lambda.min, 3), row.names = NULL)
  }
  gl <- do.call(rbind, rows); save_csv(gl, "07_glmnet_terms.csv")
  show(gl, 60)
  assign("gl", gl, envir = .GlobalEnv)
})

# ---- bootstrap metrics --------------------------------------------------------
main <- list(); delta <- list(); sg <- list()
section("3. Metrics with bootstrap intervals", {
  for (l in CFG$landmarks) {
    dl <- d0[d0$landmark == l, ]; n <- nrow(dl); Ps <- Pall[[l]]; mn <- names(Ps)
    bt <- array(NA_real_, c(B, length(mn), length(HORIZONS), 3),
                dimnames = list(NULL, mn, as.character(HORIZONS), c("AUC", "sBrier", "OE")))
    sg_levels <- list(); j2 <- match(SG_H, HORIZONS)
    if (l %in% SG_LM) for (v in SG_VARS) {
      lv <- as.character(dl[[v]])
      for (u in setdiff(unique(lv), c(NA, "Unknown", "missing"))) {
        s <- which(lv == u)
        if (length(s) >= MIN_N_SUBGROUP && sum(dl$ev[s] == 1 & dl$time[s] <= SG_H) >= MIN_EV_SUBGROUP)
          sg_levels[[length(sg_levels) + 1]] <- list(var = v, level = u, lv = lv)
      }
    }
    sg_bt <- lapply(mn, function(x) matrix(NA_real_, B, length(sg_levels))); names(sg_bt) <- mn
    set.seed(CFG$seed + 200 + match(l, CFG$landmarks))
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      for (mname in mn) for (j in seq_along(HORIZONS))
        bt[b, mname, j, ] <- met3(dl$time[idx], dl$ev[idx], Ps[[mname]][idx, j], HORIZONS[j])
      if (length(sg_levels)) for (mname in mn) for (q in seq_along(sg_levels)) {
        u <- sg_levels[[q]]; s <- idx[which(u$lv[idx] == u$level)]
        if (length(s) >= 5)
          sg_bt[[mname]][b, q] <- km_risk(dl$time[s], dl$ev[s], SG_H) / mean(Ps[[mname]][s, j2], na.rm = TRUE)
      }
      if (b %% 100 == 0) cat("  ", l, ": bootstrap", b, "of", B, "\n")
    }
    for (mname in mn) for (j in seq_along(HORIZONS)) {
      pt <- met3(dl$time, dl$ev, Ps[[mname]][, j], HORIZONS[j])
      ci <- lapply(c("AUC", "sBrier", "OE"), function(k) q95(bt[, mname, j, k]))
      main[[length(main) + 1]] <- data.frame(landmark = l, model = mname, horizon = HORIZONS[j], n = n,
        events_by_h = sum(dl$ev == 1 & dl$time <= HORIZONS[j]),
        AUC = pt[["AUC"]], AUC_lcl = ci[[1]][1], AUC_ucl = ci[[1]][2],
        scaled_Brier = pt[["sBrier"]], sBrier_lcl = ci[[2]][1], sBrier_ucl = ci[[2]][2],
        O_over_E = pt[["OE"]], OE_lcl = ci[[3]][1], OE_ucl = ci[[3]][2])
      if (mname != REF) {
        pr <- met3(dl$time, dl$ev, Ps[[REF]][, j], HORIZONS[j])[["AUC"]]
        cd <- q95(bt[, mname, j, "AUC"] - bt[, REF, j, "AUC"])
        delta[[length(delta) + 1]] <- data.frame(landmark = l, model = mname, horizon = HORIZONS[j],
          delta_AUC_vs_ref = pt[["AUC"]] - pr, lcl = cd[1], ucl = cd[2])
      }
    }
    for (mname in mn) for (q in seq_along(sg_levels)) {
      u <- sg_levels[[q]]; s <- which(u$lv == u$level); ci <- q95(sg_bt[[mname]][, q])
      ob <- km_risk(dl$time[s], dl$ev[s], SG_H); ex <- mean(Ps[[mname]][s, j2], na.rm = TRUE)
      sg[[length(sg) + 1]] <- data.frame(landmark = l, model = mname, horizon = SG_H, variable = u$var, level = u$level,
        n = length(s), events_by_h = sum(dl$ev[s] == 1 & dl$time[s] <= SG_H),
        km_observed = ob, mean_pred = ex, O_over_E = ob / ex, lcl = ci[1], ucl = ci[2])
    }
  }
  rd <- function(x) { nn <- vapply(x, is.numeric, logical(1)); x[nn] <- lapply(x[nn], round, 3); x }
  main_df <- rd(do.call(rbind, main)); delta_df <- rd(do.call(rbind, delta)); sg_df <- rd(do.call(rbind, sg))
  main_df$model <- factor(main_df$model, levels = sapply(MODELS, `[[`, "name"))
  main_df <- main_df[order(match(main_df$landmark, CFG$landmarks), main_df$model, main_df$horizon), ]
  main_df$model <- as.character(main_df$model)
  save_csv(main_df, "07_model_metrics.csv"); save_csv(delta_df, "07_auc_gain_vs_ref.csv"); save_csv(sg_df, "07_subgroup_oe_ci.csv")
  assign("main_df", main_df, envir = .GlobalEnv)
  cat("\nTwo-year horizon (95% bootstrap intervals; all horizons in 07_model_metrics.csv):\n")
  show(main_df[main_df$horizon == 2, c("landmark", "model", "events_by_h", "AUC", "AUC_lcl", "AUC_ucl",
                                       "scaled_Brier", "O_over_E", "OE_lcl", "OE_ucl")], 60)
  cat("\nAUC gain over", REF, "(paired bootstrap), 2-year horizon:\n")
  show(delta_df[delta_df$horizon == 2, ], 40)
  cat("\nSubgroup O/E at baseline, 2-year horizon (protocol, race, ethnicity):\n")
  show(sg_df[sg_df$landmark == "bl" & sg_df$variable %in% c("origprot", "race_bw", "hisp"),
             c("model", "variable", "level", "n", "events_by_h", "O_over_E", "lcl", "ucl")], 90)
})

# ---- save predictions and checks ------------------------------------------------
section("4. Save predictions and run checks", {
  saveRDS(list(d = d0[, c("RID", "landmark", "t0", "time", "ev", "fold", "sex", "race_bw", "hisp", "edu_grp",
                          "apoe4_grp", "age_grp", "origprot", "era")], P = Pall, horizons = HORIZONS),
          file.path(CFG$dir_out, "interim", "oof_predictions.rds"))
  cat("Saved output/interim/oof_predictions.rds\n\n")

  ref5 <- tryCatch(readr::read_csv(file.path(CFG$dir_out, "audit", "05_tier_metrics.csv"), show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(ref5)) {
    a <- main_df[main_df$model == "Cox-B", c("landmark", "horizon", "AUC")]
    b5 <- as.data.frame(ref5[ref5$tier == "B", c("landmark", "horizon", "AUC")])
    m <- merge(a, b5, by = c("landmark", "horizon"), suffixes = c("_07", "_05"))
    check("Cox-B reproduces the 05 Tier B AUC at every landmark and horizon (within 0.002)",
          nrow(m) == 12 && all(abs(m$AUC_07 - m$AUC_05) < 0.002), sprintf("(largest difference %.4f)", max(abs(m$AUC_07 - m$AUC_05))))
  } else check("05 results found for the regression check", FALSE)
  allP <- unlist(lapply(Pall, function(x) lapply(x, identity)), recursive = FALSE)
  check("all predicted risks lie between 0 and 1", all(sapply(allP, function(P) all(P >= 0 & P <= 1, na.rm = TRUE))))
  check("risk never decreases with the horizon, for every model", all(sapply(allP, function(P)
    all(P[, 1] <= P[, 2] + 1e-9 & P[, 2] <= P[, 3] + 1e-9, na.rm = TRUE))))
  check("fewer than 1% of out-of-fold predictions are missing, for every model",
        all(sapply(allP, function(P) mean(is.na(P[, 1])) < 0.01)),
        sprintf("(largest share %.2f%%)", 100 * max(sapply(allP, function(P) mean(is.na(P[, 1]))))))
  check("elastic net keeps at least one term at every landmark", all(tapply(gl$kept, gl$landmark, any)))
  set.seed(CFG$seed + 1)
  dl <- d0[d0$landmark == "bl", ]; idx <- sample.int(nrow(dl)); dl$time <- dl$time[idx]; dl$ev <- dl$ev[idx]
  for (nm in c("GLMNET-C", "RSF-C")) {
    mm <- MODELS[[which(sapply(MODELS, `[[`, "name") == nm)]]
    P <- predict_model(dl, mm); au <- met3(dl$time, dl$ev, P[, 2], 2)[["AUC"]]
    check(paste0("shuffled outcomes give AUC near 0.5 for ", nm, " at bl"), abs(au - 0.5) < 0.07, sprintf("(AUC = %.3f)", au))
  }
  cat(sprintf("\n%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  save_csv(tests, "07_tests.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
