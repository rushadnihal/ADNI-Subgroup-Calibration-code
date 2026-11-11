# utils_eval.R -------------------------------------------------------------
# Shared evaluation helpers for scripts 06 onward (04 and 05 hold identical
# copies of cens_G, km_risk, met3 and run_cv; they are left untouched).
# Sourced by the numbered scripts; defines functions only. Needs CFG (00_setup.R).
# -------------------------------------------------------------------------
if (!exists("K")) K <- 5
if (!exists("HORIZONS")) HORIZONS <- CFG$horizons_years

# Kaplan-Meier censoring distribution G(t) as a function (floored at 0.05).
cens_G <- function(time, ev) {
  km <- survival::survfit(survival::Surv(time, 1 - ev) ~ 1)
  if (length(km$time) == 0) return(function(t) rep(1, length(t)))
  sf <- stepfun(km$time, c(1, km$surv))
  function(t) pmax(sf(t), 0.05)
}
# Kaplan-Meier risk of the event by horizon h.
km_risk <- function(time, ev, h) {
  if (sum(ev) == 0) return(0)
  1 - summary(survival::survfit(survival::Surv(time, ev) ~ 1), times = h, extend = TRUE)$surv
}
# IPCW (Hajek) risk of the event by h. GT = G(T_i | x_i), Gh = G(h | x_i); idx selects a subgroup.
ipcw_risk <- function(time, ev, GT, Gh, h, idx = rep(TRUE, length(time))) {
  case <- ev == 1 & time <= h; ctrl <- time > h
  w <- ifelse(case, 1 / GT, ifelse(ctrl, 1 / Gh, 0))
  den <- sum(w[(case | ctrl) & idx])
  if (den == 0) return(NA_real_)
  sum(w[case & idx]) / den
}
# AUC, scaled Brier and O/E at horizon h from predicted risks p (IPCW, KM censoring).
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
# Out-of-fold Cox predictions using the fold column in d (subject-level folds).
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

# Subgroup columns used across scripts.
add_groups <- function(df) {
  df$race_bw <- ifelse(df$race %in% c("Black", "White"), df$race, NA_character_)
  df$era <- factor(ifelse(df$origprot == "ADNI1", "ADNI1",
                   ifelse(df$origprot %in% c("ADNIGO", "ADNI2"), "ADNIGO/2", "ADNI3/4")),
                   levels = c("ADNI1", "ADNIGO/2", "ADNI3/4"))
  df
}
# Analysis table: outcomes rows flagged by in_col, joined to visit predictors and folds.
make_d0 <- function(oc, vis, folds, in_col = "in_set") {
  vv <- vis[, c("RID", "EXAMDATE", "cdrsb", "apoe4", "edu", "adas13", "ravlt_immediate", "faq",
                "hippocampus", "entorhinal", "mid_temporal", "ventricles", "icv")]
  d <- oc[oc[[in_col]], ]
  d <- merge(d, vv, by.x = c("RID", "t0"), by.y = c("RID", "EXAMDATE"), all.x = TRUE, sort = FALSE)
  d <- merge(d, folds[, c("RID", "fold")], by = "RID", all.x = TRUE, sort = FALSE)
  d$time <- d$time_cens; d$ev <- d$status_cens; d$age <- d$age_lm
  d$male <- as.integer(d$sex == "Male")
  d$hippo <- d$hippocampus / d$icv * 1000; d$midtemp <- d$mid_temporal / d$icv * 1000
  d$entorhinal <- d$entorhinal / d$icv * 1000; d$ventricles <- d$ventricles / d$icv * 1000
  d <- add_groups(d)
  d <- d[order(match(d$landmark, CFG$landmarks), d$RID), ]
  rownames(d) <- NULL
  d
}
