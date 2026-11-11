# 14_temporal_holdout.R ------------------------------------------------------
# Two checks that back the paper's design claims.
#
# A. TEMPORAL HOLD-OUT. The era finding (ADNI1 under-predicted, ADNI3 over-predicted) comes from
#    out-of-fold predictions in which every era is in the training data. Here Cox-B (Tier B) is
#    trained on ADNI1 subjects only and applied to participants first enrolled in ADNI-GO, ADNI2
#    and ADNI3 (ORIGPROT, so ADNI1 rollovers are never in the test set; subjects are checked to be
#    disjoint). ADNI4 is excluded: its 2-year status is mostly unobservable. For each landmark
#    (bl, m12, m24), horizon and test group the script reports O/E and AUC of
#       A  the ADNI1-trained model (for ADNI1 itself: 5-fold CV within ADNI1)
#       B  the pooled out-of-fold Cox-B predictions saved by 07 (the paper's main analysis)
#    and the paired ratio A/B, with subject-level bootstrap intervals (one resample shared by A and
#    B, trained models are not refit, so training uncertainty is not included).
#    Reading: if A is far below 1 in ADNI3 as well, the over-prediction is a genuine shift away from
#    the ADNI1 development population and not an artefact of mixing eras in the folds.
#
# B. RECORD-LEVEL VERSUS SUBJECT-LEVEL SPLITS. The per-landmark analysis has one row per subject, so
#    splitting by record equals splitting by subject there. Leakage arises when landmarks are stacked
#    (a "super-model"): the same subject's rows then sit in both training and test folds. Cox-B plus a
#    landmark term is fitted on the stacked rows with (i) the subject-level folds and (ii) random
#    record-level folds, and the AUC is compared pooled and by landmark (bootstrap over subjects).
#
# Death is censoring throughout. Counts and metrics only in the CSVs.
# Writes: output/audit/14_*.csv, 14_temporal_holdout_report.txt, output/figures/14_holdout_and_split.png (+pdf)
# Run:    source("R/14_temporal_holdout.R")   (after 07; about 15 to 25 minutes)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
if (!requireNamespace("survival", quietly = TRUE))
  stop("Package 'survival' is missing. It ships with R; run install.packages('survival').")
K <- 5; HORIZONS <- CFG$horizons_years
source(if (file.exists("R/utils_eval.R")) "R/utils_eval.R" else "utils_eval.R")

report <- file.path(CFG$dir_out, "audit", "14_temporal_holdout_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Temporal hold-out and split check | run:", format(Sys.time()), "\n")

# ---- settings -------------------------------------------------------------
B            <- 1000
B2           <- 500
MIN_EV       <- 5
MIN_TRAIN_EV <- 30
FB      <- c("age", "male", "edu", "apoe4", "cdrsb", "adas13", "ravlt_immediate", "faq")
LM_H    <- c("bl", "m12", "m24")
LM_AUC  <- c("bl", "m12")
HS_AUC  <- c(2, 3)
GRP     <- c("ADNI1", "ADNIGO/2", "ADNI3", "ADNI2/GO/3")
NEWP    <- c("ADNIGO", "ADNI2", "ADNI3")
PAL     <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
# ---------------------------------------------------------------------------

tests <- data.frame(test = character(), result = character(), detail = character(), stringsAsFactors = FALSE)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s %s\n", if (ok) "PASS" else "FAIL", name, detail))
  tests <<- rbind(tests, data.frame(test = name, result = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE))
  invisible(ok)
}
rd <- function(x, k = 3) { nn <- vapply(x, is.numeric, logical(1)); x[nn] <- lapply(x[nn], round, k); x }

km_fast <- function(time, ev, h) {
  o <- order(time); t <- time[o]; e <- ev[o]
  te <- t[e == 1 & t <= h]
  if (!length(te)) return(0)
  r <- rle(te)
  nr <- length(t) - findInterval(r$values, t, left.open = TRUE)
  1 - prod(1 - r$lengths / nr)
}
masks_of <- function(op) {
  m <- list(op == "ADNI1", op %in% c("ADNIGO", "ADNI2"), op == "ADNI3", op %in% NEWP)
  names(m) <- GRP
  m
}
oe_stat <- function(tm, ev, P, mask, j, h) {
  idx <- which(mask & !is.na(P[, j]))
  if (length(idx) < 2) return(NA_real_)
  km_fast(tm[idx], ev[idx], h) / mean(P[idx, j])
}
auc_stat <- function(tm, ev, P, mask, j, h) {
  idx <- which(mask & !is.na(P[, j]))
  if (length(idx) < 20 || sum(ev[idx] == 1 & tm[idx] <= h) < 2) return(NA_real_)
  unname(met3(tm[idx], ev[idx], P[idx, j], h)["AUC"])
}
# Cox fitted on dtr, risks at HORIZONS for dte (median imputation from the training rows).
hold_pred <- function(dtr, dte, feats) {
  med <- sapply(feats, function(f) median(dtr[[f]], na.rm = TRUE))
  imp <- function(x) { for (f in feats) x[[f]][is.na(x[[f]])] <- med[[f]]; x }
  dtr <- imp(dtr); dte <- imp(dte)
  fml <- stats::as.formula(paste("survival::Surv(time, ev) ~", paste(feats, collapse = " + ")))
  fit <- survival::coxph(fml, data = dtr)
  sf  <- survival::survfit(fit, newdata = dte)
  P <- sapply(HORIZONS, function(h) 1 - as.numeric(summary(sf, times = h, extend = TRUE)$surv))
  list(P = matrix(P, nrow = nrow(dte)), fit = fit)
}

oc    <- readRDS(file.path(CFG$dir_out, "interim", "outcomes.rds"))
vis   <- readRDS(file.path(CFG$dir_out, "interim", "visits.rds"))
folds <- readRDS(file.path(CFG$dir_out, "interim", "folds.rds"))
o07   <- readRDS(file.path(CFG$dir_out, "interim", "oof_predictions.rds"))

d0 <- NULL; PA <- list(); PB <- list(); hr_df <- NULL; train_info <- NULL
oe_df <- ratio_df <- auc_df <- NULL
d_all <- NULL; leak_df <- NULL; diff_df <- NULL; share_s <- NA_real_; share_r <- NA_real_

section("0. Data and the pooled out-of-fold predictions from 07", {
  d0 <- make_d0(oc, vis, folds, "in_set")
  for (l in CFG$landmarks) PB[[l]] <- o07$P[[l]][["Cox-B"]]
  cat("Modelling rows:", nrow(d0), "| subjects:", length(unique(d0$RID)), "\n")
  print(table(d0$landmark, d0$origprot))
})

section("A1. Train on ADNI1, predict ADNI-GO, ADNI2 and ADNI3", {
  rows <- list(); hrs <- list()
  for (l in LM_H) {
    dl <- d0[d0$landmark == l, ]; n <- nrow(dl)
    tr <- dl$origprot == "ADNI1"; te <- dl$origprot %in% NEWP
    P <- matrix(NA_real_, n, length(HORIZONS))
    ev_tr <- sum(dl$ev[tr])
    rows[[length(rows) + 1]] <- data.frame(landmark = l, n_train = sum(tr), events_train = ev_tr, n_test = sum(te),
                                           n_adni4_not_used = sum(dl$origprot == "ADNI4"),
                                           subjects_shared = length(intersect(dl$RID[tr], dl$RID[te])))
    if (ev_tr >= MIN_TRAIN_EV) {
      P[tr, ] <- run_cv(dl[tr, ], FB)$P
      hp <- hold_pred(dl[tr, ], dl[te, ], FB)
      P[te, ] <- hp$P
      f <- hp$fit; ci <- suppressWarnings(summary(f)$conf.int)
      hrs[[length(hrs) + 1]] <- data.frame(landmark = l, term = rownames(ci), HR = ci[, "exp(coef)"],
                                           lcl = ci[, "lower .95"], ucl = ci[, "upper .95"], row.names = NULL)
    }
    PA[[l]] <- P
  }
  train_info <- do.call(rbind, rows)
  hr_df <- if (length(hrs)) do.call(rbind, hrs) else NULL
  save_csv(train_info, "14_train_test_sizes.csv"); show(train_info)
  if (!is.null(hr_df)) { save_csv(rd(hr_df), "14_adni1_trained_hr.csv"); cat("\nHazard ratios of the ADNI1-trained model saved to 14_adni1_trained_hr.csv\n") }
})

section("A2. O/E and AUC, ADNI1-trained versus pooled out-of-fold (this is the slow step)", {
  set.seed(CFG$seed + 1700)
  nG <- length(GRP); nH <- length(HORIZONS); nA <- length(HS_AUC)
  SRC <- c("A", "B")
  oe_rows <- list(); ra_rows <- list(); au_rows <- list()
  for (l in LM_H) {
    dl <- d0[d0$landmark == l, ]; n <- nrow(dl)
    Ps <- list(PA[[l]], PB[[l]])
    mk <- masks_of(dl$origprot)
    pt <- array(NA_real_, c(nG, nH, 2)); pn <- array(NA_integer_, c(nG, nH, 2)); pe <- array(NA_integer_, c(nG, nH, 2))
    kmr <- array(NA_real_, c(nG, nH, 2)); mp <- array(NA_real_, c(nG, nH, 2))
    for (g in 1:nG) for (j in 1:nH) for (s in 1:2) {
      idx <- which(mk[[g]] & !is.na(Ps[[s]][, j]))
      pn[g, j, s] <- length(idx)
      if (length(idx) >= 2) {
        h <- HORIZONS[j]
        kmr[g, j, s] <- km_fast(dl$time[idx], dl$ev[idx], h); mp[g, j, s] <- mean(Ps[[s]][idx, j])
        pt[g, j, s] <- kmr[g, j, s] / mp[g, j, s]
        pe[g, j, s] <- sum(dl$ev[idx] == 1 & dl$time[idx] <= h)
      }
    }
    doa <- l %in% LM_AUC
    pa <- array(NA_real_, c(nG, nA, 2))
    if (doa) for (g in 1:nG) for (jj in 1:nA) for (s in 1:2)
      pa[g, jj, s] <- auc_stat(dl$time, dl$ev, Ps[[s]], mk[[g]], match(HS_AUC[jj], HORIZONS), HS_AUC[jj])
    bs  <- array(NA_real_, c(B, nG, nH, 2))
    bsa <- array(NA_real_, c(B, nG, nA, 2))
    for (b in 1:B) {
      ii <- sample.int(n, replace = TRUE)
      tm <- dl$time[ii]; ev <- dl$ev[ii]; mki <- masks_of(dl$origprot[ii])
      Pi <- list(PA[[l]][ii, , drop = FALSE], PB[[l]][ii, , drop = FALSE])
      for (g in 1:nG) for (s in 1:2) {
        for (j in 1:nH) bs[b, g, j, s] <- oe_stat(tm, ev, Pi[[s]], mki[[g]], j, HORIZONS[j])
        if (doa) for (jj in 1:nA) bsa[b, g, jj, s] <- auc_stat(tm, ev, Pi[[s]], mki[[g]], match(HS_AUC[jj], HORIZONS), HS_AUC[jj])
      }
    }
    for (g in 1:nG) for (j in 1:nH) for (s in 1:2) {
      ci <- q95(bs[, g, j, s])
      oe_rows[[length(oe_rows) + 1]] <- data.frame(landmark = l, group = GRP[g], horizon = HORIZONS[j],
        source = if (s == 1) (if (g == 1) "A within-ADNI1 5-fold CV" else "A trained on ADNI1") else "B pooled out-of-fold (07)",
        n = pn[g, j, s], events_by_h = pe[g, j, s], km_risk = kmr[g, j, s], mean_pred = mp[g, j, s],
        O_over_E = pt[g, j, s], lcl = ci[1], ucl = ci[2], low_events = !is.na(pe[g, j, s]) && pe[g, j, s] < MIN_EV,
        stringsAsFactors = FALSE)
    }
    for (g in 2:nG) for (j in 1:nH) {
      r <- bs[, g, j, 1] / bs[, g, j, 2]; ci <- q95(r)
      ra_rows[[length(ra_rows) + 1]] <- data.frame(landmark = l, group = GRP[g], horizon = HORIZONS[j],
        events_by_h = pe[g, j, 1], OE_trained_on_ADNI1 = pt[g, j, 1], OE_pooled = pt[g, j, 2],
        ratio_A_over_B = pt[g, j, 1] / pt[g, j, 2], lcl = ci[1], ucl = ci[2], stringsAsFactors = FALSE)
    }
    if (doa) for (g in 1:nG) for (jj in 1:nA) {
      for (s in 1:2) {
        ci <- q95(bsa[, g, jj, s])
        au_rows[[length(au_rows) + 1]] <- data.frame(landmark = l, group = GRP[g], horizon = HS_AUC[jj],
          source = if (s == 1) (if (g == 1) "A within-ADNI1 5-fold CV" else "A trained on ADNI1") else "B pooled out-of-fold (07)",
          events_by_h = pe[g, match(HS_AUC[jj], HORIZONS), s], AUC = pa[g, jj, s], lcl = ci[1], ucl = ci[2], stringsAsFactors = FALSE)
      }
    }
    cat("  landmark", l, "done\n")
  }
  oe_df <- do.call(rbind, oe_rows); ratio_df <- do.call(rbind, ra_rows)
  auc_df <- if (length(au_rows)) do.call(rbind, au_rows) else NULL
  save_csv(rd(oe_df), "14_holdout_oe.csv"); save_csv(rd(ratio_df), "14_holdout_oe_ratio.csv")
  if (!is.null(auc_df)) save_csv(rd(auc_df), "14_holdout_auc.csv")
  cat("\nO/E at baseline, 2-year horizon (A = ADNI1-trained, B = pooled out-of-fold):\n")
  show(rd(oe_df[oe_df$landmark == "bl" & oe_df$horizon == 2, c("group", "source", "n", "events_by_h", "O_over_E", "lcl", "ucl")]), 20)
  cat("\nPaired ratio A / B (1 means the ADNI1-trained model is calibrated like the pooled one):\n")
  show(rd(ratio_df[ratio_df$horizon == 2, ]), 20)
  if (!is.null(auc_df)) { cat("\nAUC:\n"); show(rd(auc_df[auc_df$horizon == 2, ]), 30) }
  cat("\nGroups with fewer than", MIN_EV, "events by the horizon are flagged low_events and should not be interpreted.\n")
})

section("B. Record-level versus subject-level splits on stacked landmark rows", {
  set.seed(CFG$seed + 1800)
  d_all <<- d0[!is.na(d0$fold), ]
  d_all$lmf <- factor(d_all$landmark, levels = CFG$landmarks)
  FBL <- c(FB, "lmf")
  n <- nrow(d_all)
  P_s <- run_cv(d_all, FBL)$P
  d_r <- d_all; d_r$fold <- sample(rep(1:K, length.out = n))
  P_r <- run_cv(d_r, FBL)$P
  share <- function(d) {
    out <- logical(nrow(d))
    for (k in 1:K) { te <- d$fold == k; out[te] <- d$RID[te] %in% d$RID[!te] }
    mean(out)
  }
  share_s <<- share(d_all); share_r <<- share(d_r)
  cat("Rows:", n, "| subjects:", length(unique(d_all$RID)), "\n")
  cat(sprintf("Share of test rows whose subject also has a row in the training folds: subject-level %.3f, record-level %.3f\n", share_s, share_r))
  scopes <- c("pooled", CFG$landmarks); H2 <- 2; H3 <- 3
  spec <- rbind(data.frame(scope = "pooled", h = c(H2, H3), stringsAsFactors = FALSE),
                data.frame(scope = CFG$landmarks, h = H2, stringsAsFactors = FALSE))
  sel <- function(lmk, scope) if (scope == "pooled") rep(TRUE, length(lmk)) else lmk == scope
  m_one <- function(tm, ev, P, lmk, scope, h) {
    ix <- which(sel(lmk, scope))
    if (length(ix) < 20 || sum(ev[ix] == 1 & tm[ix] <= h) < 2) return(c(AUC = NA_real_, sBrier = NA_real_, OE = NA_real_))
    met3(tm[ix], ev[ix], P[ix, match(h, HORIZONS)], h)
  }
  PP <- list(subject = P_s, record = P_r)
  pt <- array(NA_real_, c(nrow(spec), 2, 3), dimnames = list(NULL, names(PP), c("AUC", "sBrier", "OE")))
  for (r in seq_len(nrow(spec))) for (v in 1:2)
    pt[r, v, ] <- m_one(d_all$time, d_all$ev, PP[[v]], d_all$landmark, spec$scope[r], spec$h[r])
  rows_by <- split(seq_len(n), d_all$RID); rids <- names(rows_by)
  bs <- array(NA_real_, c(B2, nrow(spec), 2))
  for (b in 1:B2) {
    ii <- unlist(rows_by[sample(rids, replace = TRUE)], use.names = FALSE)
    tm <- d_all$time[ii]; ev <- d_all$ev[ii]; lm <- d_all$landmark[ii]
    for (r in seq_len(nrow(spec))) for (v in 1:2)
      bs[b, r, v] <- m_one(tm, ev, PP[[v]][ii, , drop = FALSE], lm, spec$scope[r], spec$h[r])["AUC"]
  }
  lk <- list(); dfl <- list()
  for (r in seq_len(nrow(spec))) {
    ix <- which(sel(d_all$landmark, spec$scope[r]))
    for (v in 1:2) {
      ci <- q95(bs[, r, v])
      lk[[length(lk) + 1]] <- data.frame(split = names(PP)[v], scope = spec$scope[r], horizon = spec$h[r], n = length(ix),
        events_by_h = sum(d_all$ev[ix] == 1 & d_all$time[ix] <= spec$h[r]), AUC = pt[r, v, "AUC"], lcl = ci[1], ucl = ci[2],
        O_over_E = pt[r, v, "OE"], stringsAsFactors = FALSE)
    }
    ci <- q95(bs[, r, 2] - bs[, r, 1])
    dfl[[length(dfl) + 1]] <- data.frame(scope = spec$scope[r], horizon = spec$h[r],
      AUC_subject = pt[r, 1, "AUC"], AUC_record = pt[r, 2, "AUC"], inflation = pt[r, 2, "AUC"] - pt[r, 1, "AUC"], lcl = ci[1], ucl = ci[2],
      stringsAsFactors = FALSE)
  }
  leak_df <<- do.call(rbind, lk); diff_df <<- do.call(rbind, dfl)
  save_csv(rd(leak_df), "14_split_auc.csv"); save_csv(rd(diff_df), "14_split_auc_inflation.csv")
  cat("\nAUC inflation from record-level splitting (record minus subject), bootstrap over subjects:\n")
  show(rd(diff_df), 20)
  cat("\nThis is a property of stacking landmarks. The per-landmark models in 04 to 13 have one row per subject, so they are not exposed.\n")
})

section("C. Figure", {
  fdir <- file.path(CFG$dir_out, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  ci_v <- function(x, l, u, lo, hi, col) {
    if (!is.finite(l) || !is.finite(u)) return(invisible())
    segments(x, max(l, lo), x, min(u, hi), col = col, lwd = 1.5)
    if (l < lo) points(x, lo, pch = "v", col = col, cex = 0.7)
    if (u > hi) points(x, hi, pch = "^", col = col, cex = 0.7)
  }
  fig <- function() {
    layout(matrix(1:3, 1), widths = c(1, 1, 1.15)); par(mar = c(5, 4.2, 3, 1))
    for (l in c("bl", "m12")) {
      plot(NA, xlim = c(0.5, length(GRP) + 0.5), ylim = c(0.2, 3), log = "y", xaxt = "n", yaxt = "n",
           xlab = "Test group", ylab = "Observed / expected (2-year horizon)", main = paste("Landmark", l))
      axis(1, at = seq_along(GRP), labels = GRP, cex.axis = 0.75)
      axis(2, at = c(0.25, 0.5, 0.67, 1, 1.5, 2, 3), labels = c("0.25", "0.5", "0.67", "1", "1.5", "2", "3"), las = 1)
      abline(h = 1, lty = 2, col = "grey40")
      for (s in 1:2) for (g in seq_along(GRP)) {
        src <- if (s == 1) (if (g == 1) "A within-ADNI1 5-fold CV" else "A trained on ADNI1") else "B pooled out-of-fold (07)"
        r <- oe_df[oe_df$landmark == l & oe_df$horizon == 2 & oe_df$group == GRP[g] & oe_df$source == src, ]
        if (nrow(r) != 1 || is.na(r$events_by_h) || r$events_by_h < MIN_EV) next
        x <- g + (s - 1.5) * 0.25
        ci_v(x, r$lcl, r$ucl, 0.2, 3, PAL[s]); points(x, r$O_over_E, pch = 16, col = PAL[s])
      }
      legend("topright", legend = c("Trained on ADNI1 (ADNI1: 5-fold CV)", "Pooled out-of-fold"), col = PAL[1:2], pch = 16, bty = "n", cex = 0.65)
    }
    par(mar = c(7, 4.8, 3, 1))
    plot(NA, xlim = c(0.5, nrow(diff_df) + 0.5), ylim = range(c(0, diff_df$lcl, diff_df$ucl, diff_df$inflation), na.rm = TRUE) * c(1, 1.1),
         xaxt = "n", xlab = "", ylab = "AUC, record-level minus subject-level split", main = "Stacked landmarks")
    title(xlab = "Scope (AUC horizon)", line = 5.8)
    lab <- paste0(diff_df$scope, " (", diff_df$horizon, "y)")
    axis(1, at = seq_len(nrow(diff_df)), labels = lab, cex.axis = 0.7, las = 2)
    abline(h = 0, lty = 2, col = "grey40")
    for (i in seq_len(nrow(diff_df))) { ci_v(i, diff_df$lcl[i], diff_df$ucl[i], -1, 1, PAL[3]); points(i, diff_df$inflation[i], pch = 16, col = PAL[3]) }
  }
  ok <- tryCatch({
    grDevices::png(file.path(fdir, "14_holdout_and_split.png"), width = 11, height = 4.8, units = "in", res = 200); fig(); grDevices::dev.off()
    grDevices::pdf(file.path(fdir, "14_holdout_and_split.pdf"), width = 11, height = 4.8); fig(); grDevices::dev.off(); TRUE
  }, error = function(e) { while (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE); cat("  figure not written:", conditionMessage(e), "\n"); FALSE })
  assign("fig_ok14", ok, envir = .GlobalEnv)
  cat("Figure:", if (ok) "14_holdout_and_split.png and .pdf" else "none", "\n")
})

section("D. Checks", {
  ok_data <- !is.null(d0) && length(PA) == length(LM_H) && !is.null(oe_df)
  check("data, predictions and tables were built", ok_data)
  if (ok_data) {
    key <- function(x) paste(x$RID, x$landmark)
    check("modelling rows line up with the rows saved by 07", identical(key(d0), key(o07$d)))
    check("no subject is in both the ADNI1 training set and a test group", all(train_info$subjects_shared == 0))
    miss <- vapply(LM_H, function(l) { dl <- d0[d0$landmark == l, ]; sum(is.na(PA[[l]][dl$origprot %in% c("ADNI1", NEWP), 1])) }, numeric(1))
    check("predictions exist for every ADNI1 and test row", all(miss == 0), sprintf("(missing: %s)", paste(miss, collapse = ", ")))
    a4 <- vapply(LM_H, function(l) { dl <- d0[d0$landmark == l, ]; all(is.na(PA[[l]][dl$origprot == "ADNI4", 1])) }, logical(1))
    check("ADNI4 rows are not used as test rows", all(a4))
    dl <- d0[d0$landmark == "bl", ]
    check("km_fast reproduces survfit Kaplan-Meier risk (within 1e-9)", abs(km_fast(dl$time, dl$ev, 2) - km_risk(dl$time, dl$ev, 2)) < 1e-9)
    ev2 <- function(g) { x <- oe_df[oe_df$landmark == "bl" & oe_df$horizon == 2 & oe_df$group == g & grepl("^B", oe_df$source), ]; x$events_by_h }
    check("baseline events by 2 years: ADNI1 118, ADNIGO/2 56, ADNI3 14", identical(c(ev2("ADNI1"), ev2("ADNIGO/2"), ev2("ADNI3")), c(118L, 56L, 14L)),
          sprintf("(%s)", paste(c(ev2("ADNI1"), ev2("ADNIGO/2"), ev2("ADNI3")), collapse = ", ")))
    b3 <- oe_df[oe_df$landmark == "bl" & oe_df$horizon == 2 & oe_df$group == "ADNI3" & grepl("^B", oe_df$source), ]
    check("pooled ADNI3 O/E at baseline, 2 years matches 08 (0.502)", nrow(b3) == 1 && abs(b3$O_over_E - 0.502) < 0.01, sprintf("(%.3f)", b3$O_over_E))
    b1 <- oe_df[oe_df$landmark == "bl" & oe_df$horizon == 2 & oe_df$group == "ADNI1" & grepl("^B", oe_df$source), ]
    check("pooled ADNI1 O/E at baseline, 2 years matches 08 (1.34)", nrow(b1) == 1 && abs(b1$O_over_E - 1.34) < 0.02, sprintf("(%.3f)", b1$O_over_E))
    use <- !oe_df$low_events & !is.na(oe_df$low_events)
    check("intervals are finite for every group with enough events", all(is.finite(oe_df$lcl[use]) & is.finite(oe_df$ucl[use])), sprintf("(%d rows)", sum(use)))
    inside <- oe_df$O_over_E >= oe_df$lcl - 0.15 & oe_df$O_over_E <= oe_df$ucl + 0.15
    check("O/E point estimates sit inside their intervals (within 0.15)", all(inside[use], na.rm = TRUE))
    check("subject-level folds never split a subject; record-level folds do", isTRUE(share_s == 0) && isTRUE(share_r > 0.3),
          sprintf("(%.3f and %.3f)", share_s, share_r))
    check("split AUC table was built", !is.null(diff_df) && all(is.finite(diff_df$inflation)))
    check("figure written", isTRUE(get0("fig_ok14", ifnotfound = FALSE)) && file.exists(file.path(CFG$dir_out, "figures", "14_holdout_and_split.png")))
  }
  cat(sprintf("\n%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  save_csv(tests, "14_tests.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
