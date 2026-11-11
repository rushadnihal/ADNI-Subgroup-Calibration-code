# 12_paper_figures.R ---------------------------------------------------------
# Paper-ready figures, redrawn from the CSVs written by scripts 08 to 11. Nothing is
# refit and nothing is recomputed except the events-needed formula used for reference
# lines. Each figure is written as PNG (200 dpi) and PDF.
#
# Changes from the first-pass figures in 08 to 11:
#   - intervals that run past the axis are clipped with an arrow mark, not cut silently
#   - decision curves show the net-benefit gain over the best default strategy with its
#     interval, so "interval above 0" can be read off the figure
#   - plasma results are shown as paired AUC differences, not overlapping AUC intervals
#   - subgroup O/E estimates are drawn over the range ADNI could not distinguish from 1
#   - events-needed figure uses dots, with all reference lines labelled
#   - the KM versus Aalen-Johansen figure is one panel with a correct axis label
#
# Reads:  output/audit/08_oe_over_time.csv, 08_oe_contrasts.csv, 09_decision_curve.csv,
#         09_conformal_by_subgroup.csv, 10_plasma_auc_gain.csv, 10_cr_risk_comparison.csv,
#         11_subgroup_events_vs_needed.csv
# Writes: output/figures/12_*.png and 12_*.pdf, output/audit/12_figures_report.txt,
#         output/audit/12_tests.csv
# Run:    source("R/12_paper_figures.R")   (after 08 to 11; under a minute)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")

report <- file.path(CFG$dir_out, "audit", "12_figures_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Paper figures | run:", format(Sys.time()), "\n")

# ---- settings -------------------------------------------------------------
LMS   <- c("bl", "m12", "m24", "m36")
M4    <- c("Cox-B", "Cox-C", "GLMNET-C", "RSF-C")
ERAS  <- c("ADNI1", "ADNI2", "ADNI3")   # ADNI-GO (3 events by 2 years) and ADNI4 (not estimable at 2 years) are left out of the era figure
VARS  <- c("all", "sex", "race_bw", "hisp", "edu_grp", "apoe4_grp", "age_grp", "origprot")
PAL   <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000")
ZA    <- qnorm(0.975); ZB <- qnorm(0.80)
E_needed <- function(m) (ZA + ZB * sqrt(m))^2 / (m - 1)^2
# ---------------------------------------------------------------------------

tests <- data.frame(test = character(), result = character(), detail = character(), stringsAsFactors = FALSE)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s %s\n", if (ok) "PASS" else "FAIL", name, detail))
  tests <<- rbind(tests, data.frame(test = name, result = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE))
  invisible(ok)
}

A <- function(f) file.path(CFG$dir_out, "audit", f)
read_in <- function(f, cols) {
  p <- A(f)
  if (!file.exists(p)) { cat("  missing input:", p, "\n"); return(NULL) }
  x <- utils::read.csv(p, stringsAsFactors = FALSE)
  miss <- setdiff(cols, names(x))
  if (length(miss)) { cat("  ", f, "lacks columns:", paste(miss, collapse = ", "), "\n"); return(NULL) }
  x
}

n_clip <- 0L
fig_ok <- list()
oe <- dc <- cf <- pg <- cr <- p11 <- NULL

# vertical interval on a log y axis, clipped to [lo, hi] with an arrow mark where clipped
ci_v <- function(x, l, u, lo, hi, col, lwd = 1.5) {
  if (!is.finite(l) || !is.finite(u)) return(invisible())
  segments(x, max(l, lo), x, min(u, hi), col = col, lwd = lwd)
  if (l < lo) { points(x, lo, pch = "v", col = col, cex = 0.7); n_clip <<- n_clip + 1L }
  if (u > hi) { points(x, hi, pch = "^", col = col, cex = 0.7); n_clip <<- n_clip + 1L }
}
# horizontal interval, same idea
ci_h <- function(y, l, u, lo, hi, col, lwd = 1.5) {
  if (!is.finite(l) || !is.finite(u)) return(invisible())
  segments(max(l, lo), y, min(u, hi), y, col = col, lwd = lwd)
  if (l < lo) { points(lo, y, pch = "<", col = col, cex = 0.7); n_clip <<- n_clip + 1L }
  if (u > hi) { points(hi, y, pch = ">", col = col, cex = 0.7); n_clip <<- n_clip + 1L }
}

make_fig <- function(name, w, h, fn) {
  fdir <- file.path(CFG$dir_out, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  ok <- tryCatch({
    grDevices::png(file.path(fdir, paste0(name, ".png")), width = w, height = h, units = "in", res = 200)
    fn(); grDevices::dev.off()
    grDevices::pdf(file.path(fdir, paste0(name, ".pdf")), width = w, height = h)
    fn(); grDevices::dev.off()
    TRUE
  }, error = function(e) {
    while (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE)
    cat("  FAILED", name, ":", conditionMessage(e), "\n"); FALSE
  })
  fig_ok[[name]] <<- ok
  cat("  ", name, if (ok) "written (png, pdf)" else "NOT written", "\n")
  ok
}

section("1. Inputs", {
  oe   <<- read_in("08_oe_over_time.csv", c("landmark", "model", "horizon", "variable", "level", "n", "events_by_h", "O_over_E", "lcl", "ucl"))
  dc   <<- read_in("09_decision_curve.csv", c("landmark", "model", "horizon", "threshold", "NB_model", "delta_NB", "lcl", "ucl"))
  cf   <<- read_in("09_conformal_by_subgroup.csv", c("landmark", "model", "horizon", "method", "variable", "level", "n_cases", "coverage_ipcw", "lcl", "ucl"))
  pg   <<- read_in("10_plasma_auc_gain.csv", c("landmark", "model", "horizon", "events_by_h", "delta_AUC_vs_CoxB", "lcl", "ucl"))
  cr   <<- read_in("10_cr_risk_comparison.csv", c("landmark", "variable", "level", "horizon", "difference_pp"))
  p11  <<- read_in("11_subgroup_events_vs_needed.csv", c("landmark", "variable", "level", "n", "events_by_h", "expected_events_known",
                    "detectable_OE_below_after_loss", "detectable_OE_above_after_loss"))
  inp <- list(oe = oe, dc = dc, cf = cf, pg = pg, cr = cr, p11 = p11)
  for (nm in names(inp)) cat(sprintf("  %-4s %s\n", nm, if (is.null(inp[[nm]])) "NOT LOADED" else paste(nrow(inp[[nm]]), "rows")))
})

# ---- Figure 1: whole-cohort O/E by model ------------------------------------
f1 <- function() {
  d <- oe[oe$variable == "all" & oe$model %in% M4, ]
  par(mfrow = c(1, 3), mar = c(4, 4.2, 3, 1))
  for (h in 1:3) {
    plot(NA, xlim = c(0.5, 4.5), ylim = c(0.5, 2), log = "y", xaxt = "n", yaxt = "n",
         xlab = "Landmark", ylab = "Observed / expected", main = sprintf("%d-year horizon", h))
    axis(1, at = 1:4, labels = LMS)
    axis(2, at = c(0.5, 0.67, 0.8, 1, 1.25, 1.5, 2), labels = c("0.5", "0.67", "0.8", "1", "1.25", "1.5", "2"), las = 1)
    abline(h = 1, lty = 2, col = "grey40")
    for (k in seq_along(M4)) {
      dd <- d[d$model == M4[k] & d$horizon == h, ]
      if (!nrow(dd)) next
      dd <- dd[order(match(dd$landmark, LMS)), ]
      x <- match(dd$landmark, LMS) + (k - 2) * 0.22
      lines(x, dd$O_over_E, col = PAL[k], lwd = 1)
      for (i in seq_len(nrow(dd))) ci_v(x[i], dd$lcl[i], dd$ucl[i], 0.5, 2, PAL[k])
      points(x, dd$O_over_E, pch = 16, col = PAL[k], cex = 0.9)
    }
    if (h == 1) legend("topleft", legend = M4, col = PAL[1:4], pch = 16, lwd = 1.5, bty = "n", cex = 0.8, ncol = 2)
  }
}

# ---- Figure 2: O/E by protocol era over landmarks ---------------------------
f2 <- function() {
  par(mfrow = c(1, 2), mar = c(5, 4.2, 3, 1))
  for (mod in c("Cox-B", "Cox-B+era")) {
    plot(NA, xlim = c(0.5, 4.5), ylim = c(0.2, 3), log = "y", xaxt = "n", yaxt = "n",
         xlab = "Landmark", ylab = "Observed / expected (2-year horizon)", main = mod)
    axis(1, at = 1:4, labels = LMS)
    axis(2, at = c(0.25, 0.5, 0.67, 1, 1.5, 2, 3), labels = c("0.25", "0.5", "0.67", "1", "1.5", "2", "3"), las = 1)
    abline(h = 1, lty = 2, col = "grey40")
    for (k in seq_along(ERAS)) {
      dd <- oe[oe$variable == "origprot" & oe$level == ERAS[k] & oe$model == mod & oe$horizon == 2, ]
      if (!nrow(dd)) next
      dd <- dd[order(match(dd$landmark, LMS)), ]
      x <- match(dd$landmark, LMS) + (k - 2) * 0.22
      text(x, 0.21, dd$events_by_h, cex = 0.5, col = PAL[k])
      keep <- dd$events_by_h >= 5               # points on fewer than 5 events are not drawn
      dd <- dd[keep, ]; x <- x[keep]
      if (!nrow(dd)) next
      lines(x, dd$O_over_E, col = PAL[k], lwd = 1)
      for (i in seq_len(nrow(dd))) ci_v(x[i], dd$lcl[i], dd$ucl[i], 0.2, 3, PAL[k])
      points(x, dd$O_over_E, pch = 16, col = PAL[k], cex = 0.9)
    }
    legend("topright", legend = ERAS, col = PAL[1:3], pch = 16, lwd = 1.5, bty = "n", cex = 0.75)
  }
}

# ---- Figure 3: subgroup O/E at baseline over the undetectable range ---------
f3 <- function() {
  d <- oe[oe$landmark == "bl" & oe$model == "Cox-B" & oe$horizon == 2 & oe$variable %in% VARS, ]
  d <- merge(d, p11[p11$landmark == "bl", c("variable", "level", "detectable_OE_below_after_loss", "detectable_OE_above_after_loss")],
             by = c("variable", "level"), all.x = TRUE, sort = FALSE)
  d <- d[order(match(d$variable, VARS)), ]
  n <- nrow(d); y <- rev(seq_len(n))
  lab <- ifelse(d$variable == "all", sprintf("all (%d ev)", d$events_by_h), sprintf("%s: %s (%d ev)", d$variable, d$level, d$events_by_h))
  lo <- 0.1; hi <- 4
  par(mar = c(4.2, 11, 3, 1))
  plot(NA, xlim = c(lo, hi), ylim = c(0.5, n + 2.5), log = "x", xaxt = "n", yaxt = "n", xlab = "Observed / expected, baseline, 2 years (Cox-B)", ylab = "",
       main = "Subgroup O/E and the range ADNI could not tell from 1", cex.main = 0.95)
  axis(1, at = c(0.1, 0.25, 0.5, 1, 2, 4), labels = c("0.1", "0.25", "0.5", "1", "2", "4"))
  axis(2, at = y, labels = lab, las = 1, cex.axis = 0.65)
  for (i in seq_len(n)) {
    b <- d$detectable_OE_below_after_loss[i]; a <- d$detectable_OE_above_after_loss[i]
    if (is.finite(b) && is.finite(a)) rect(max(b, lo), y[i] - 0.38, min(a, hi), y[i] + 0.38, col = "grey88", border = NA)
  }
  abline(v = 1, lty = 2, col = "grey40")
  for (i in seq_len(n)) {
    ci_h(y[i], d$lcl[i], d$ucl[i], lo, hi, PAL[1])
    if (is.finite(d$O_over_E[i])) points(d$O_over_E[i], y[i], pch = 16, col = PAL[1], cex = 0.9)
  }
  legend("topright", legend = c("Estimate and 95% interval", "O/E values with under 80% power to differ from 1"),
         pch = c(16, 15), col = c(PAL[1], "grey80"), bty = "o", bg = "white", cex = 0.6)
}

# ---- Figure 4: events available versus needed -------------------------------
f4 <- function() {
  d <- p11[p11$landmark == "bl", ]
  d <- d[order(d$expected_events_known), ]
  n <- nrow(d); y <- seq_len(n)
  lab <- ifelse(d$variable == "all", "all", paste0(d$variable, ": ", d$level))
  par(mar = c(4.2, 11, 3.5, 1))
  plot(NA, xlim = c(1, 300), ylim = c(0.5, n + 0.5), log = "x", yaxt = "n", xaxt = "n",
       xlab = "Events (log scale)", ylab = "", main = "Baseline, 2 years: events available versus events needed")
  axis(1, at = c(1, 2, 5, 10, 20, 50, 100, 200), labels = c("1", "2", "5", "10", "20", "50", "100", "200"))
  axis(2, at = y, labels = lab, las = 1, cex.axis = 0.65)
  nd <- c(0.67, 0.8, 1.25)
  for (m in nd) abline(v = E_needed(m), lty = c(3, 2, 2)[match(m, nd)], col = "grey40")
  mtext(sprintf("Expected events needed for 80%% power: O/E 0.67 = %d, 1.25 = %d, 0.8 = %d (dashed lines, left to right)",
                ceiling(E_needed(0.67)), ceiling(E_needed(1.25)), ceiling(E_needed(0.8))), side = 3, line = 0.3, cex = 0.6)
  segments(d$expected_events_known, y - 0.25, d$expected_events_known, y + 0.25, col = PAL[1], lwd = 2)
  points(d$expected_events_known, y, pch = 16, col = PAL[1], cex = 1.1)
  points(pmax(d$events_by_h, 1), y, pch = 4, col = PAL[2], lwd = 1.5)
  legend("bottomright", legend = c("Expected events, subjects with known status", "Observed events by 2 years"), pch = c(16, 4),
         col = c(PAL[1], PAL[2]), bty = "o", bg = "white", cex = 0.6)
}

# ---- Figure 5: decision curves, gain over best default ----------------------
f5 <- function() {
  par(mfrow = c(1, 4), mar = c(4, 4.2, 3, 1))
  d <- dc[dc$horizon == 2 & dc$model %in% c("Cox-B", "GLMNET-C", "RSF-C"), ]
  yl <- range(c(0, d$delta_NB, d$lcl[d$model == "Cox-B"], d$ucl[d$model == "Cox-B"]), na.rm = TRUE)
  for (l in LMS) {
    plot(NA, xlim = c(0.05, 0.5), ylim = yl, xlab = "Risk threshold", ylab = "Net benefit above best default",
         main = sprintf("%s, 2-year", l))
    abline(h = 0, col = "grey40")
    b <- d[d$landmark == l & d$model == "Cox-B", ]; b <- b[order(b$threshold), ]
    ok <- is.finite(b$lcl) & is.finite(b$ucl)
    if (any(ok)) polygon(c(b$threshold[ok], rev(b$threshold[ok])), c(b$lcl[ok], rev(b$ucl[ok])), col = "#0072B220", border = NA)
    for (k in 1:3) {
      mm <- c("Cox-B", "GLMNET-C", "RSF-C")[k]
      dd <- d[d$landmark == l & d$model == mm, ]; dd <- dd[order(dd$threshold), ]
      lines(dd$threshold, dd$delta_NB, col = PAL[k], lwd = if (k == 1) 2 else 1.2)
      points(dd$threshold, dd$delta_NB, col = PAL[k], pch = 16, cex = 0.5)
    }
    if (l == "bl") legend("topright", legend = c("Cox-B (shaded: 95% interval)", "GLMNET-C", "RSF-C"), col = PAL[1:3], lwd = c(2, 1.2, 1.2), bty = "n", cex = 0.6)
  }
}

# ---- Figure 6 and S1: conformal coverage by subgroup ------------------------
cov_fig <- function(method, title) {
  function() {
    d <- cf[cf$landmark == "bl" & cf$model == "Cox-B" & cf$horizon == 2 & cf$method == method & cf$variable %in% VARS, ]
    d <- d[order(match(d$variable, VARS)), ]
    n <- nrow(d); y <- rev(seq_len(n))
    lab <- ifelse(d$variable == "all", sprintf("all (%d ev)", d$n_cases), sprintf("%s: %s (%d ev)", d$variable, d$level, d$n_cases))
    par(mar = c(4.2, 11, 3, 1))
    plot(NA, xlim = c(0.6, 1), ylim = c(0.5, n + 0.5), yaxt = "n", xlab = "Coverage of the true outcome class (censoring-weighted)", ylab = "", main = title, cex.main = 0.85)
    axis(2, at = y, labels = lab, las = 1, cex.axis = 0.65)
    abline(v = 0.9, lty = 2, col = "grey40")
    for (i in seq_len(n)) {
      ci_h(y[i], d$lcl[i], d$ucl[i], 0.6, 1, PAL[1])
      points(d$coverage_ipcw[i], y[i], pch = if (d$n_cases[i] < 30) 1 else 16, col = PAL[1], cex = 1)
    }
    legend("topleft", legend = c("30 or more events", "Under 30 events"), pch = c(16, 1), col = PAL[1], bty = "o", bg = "white", cex = 0.6)
  }
}

# ---- Figure 7: plasma, paired AUC differences -------------------------------
f7 <- function() {
  MIN_EV_PLOT <- 10
  g <- pg[pg$events_by_h >= MIN_EV_PLOT, ]
  xr <- range(c(0, g$lcl, g$ucl), na.rm = TRUE); xr <- xr + c(-0.03, 0.03)
  par(mfrow = c(2, 3), mar = c(4.8, 10, 3, 1))
  for (l in c("bl", "m12")) for (h in 1:3) {
    dd <- pg[pg$landmark == l & pg$horizon == h, ]
    ne <- if (nrow(dd)) dd$events_by_h[1] else 0L
    dd <- dd[dd$events_by_h >= MIN_EV_PLOT, ]
    n <- nrow(dd); y <- rev(seq_len(max(n, 1)))
    plot(NA, xlim = xr, ylim = c(0.5, max(n, 1) + 0.5), yaxt = "n", xlab = "AUC difference from Cox-B on the same rows", ylab = "",
         main = sprintf("%s, %d-year (%d events)", l, h, ne), cex.main = 0.9)
    if (n == 0) { text(mean(xr), 1, sprintf("fewer than %d events: not shown", MIN_EV_PLOT), cex = 0.8); next }
    axis(2, at = y, labels = dd$model, las = 1, cex.axis = 0.6)
    abline(v = 0, lty = 2, col = "grey40")
    for (i in seq_len(n)) {
      ci_h(y[i], dd$lcl[i], dd$ucl[i], xr[1], xr[2], PAL[1])
      points(dd$delta_AUC_vs_CoxB[i], y[i], pch = 16, col = PAL[1], cex = 0.9)
    }
  }
}

# ---- Figure S2: death as a competing event ----------------------------------
fs2 <- function() {
  d <- cr[cr$variable == "all", ]
  par(mar = c(4, 4.5, 3, 1))
  plot(NA, xlim = c(0.7, 4.3), ylim = c(0, 1), xaxt = "n", xlab = "Landmark", ylab = "1 - KM minus Aalen-Johansen risk (percentage points)",
       main = "Treating death as censoring overstates risk very little")
  axis(1, at = 1:4, labels = LMS)
  abline(h = 1, lty = 3, col = "grey50"); text(0.72, 0.97, "1 percentage point", adj = c(0, 1), cex = 0.65, col = "grey40")
  for (h in 1:3) {
    dd <- d[d$horizon == h, ]; dd <- dd[order(match(dd$landmark, LMS)), ]
    lines(match(dd$landmark, LMS), dd$difference_pp, col = PAL[h], lwd = 1.8, type = "b", pch = 16)
  }
  legend("topright", legend = c("1 year", "2 years", "3 years"), col = PAL[1:3], lwd = 1.8, pch = 16, bty = "n", cex = 0.8)
}

section("2. Figures", {
  n_clip <<- 0L
  if (!is.null(oe))  { make_fig("12_fig1_overall_calibration", 10, 3.8, f1); make_fig("12_fig2_calibration_by_era", 9, 4.6, f2)
                       if (!is.null(p11)) make_fig("12_fig3_subgroup_calibration_power", 8, 6.2, f3) }
  if (!is.null(p11)) make_fig("12_fig4_events_vs_needed", 8, 6.2, f4)
  if (!is.null(dc))  make_fig("12_fig5_decision_curves", 12, 3.8, f5)
  if (!is.null(cf))  { make_fig("12_fig6_conformal_class_conditional", 8, 6.2, cov_fig("class-cond-IPCW", "Class-conditional conformal coverage by subgroup (baseline, 2 years)"))
                       make_fig("12_figS1_conformal_marginal", 8, 6.2, cov_fig("marginal-IPCW", "Marginal conformal coverage by subgroup (baseline, 2 years)")) }
  if (!is.null(pg))  make_fig("12_fig7_plasma_auc_difference", 11, 6.4, f7)
  if (!is.null(cr))  make_fig("12_figS2_competing_risk", 6, 4.5, fs2)
  cat("\nInterval ends clipped to the axis and marked with an arrow:", n_clip, "\n")
})

section("3. Suggested captions", {
  cat("Fig 1. Whole-cohort observed/expected (O/E) ratio by model, landmark and horizon, from out-of-fold predictions, with 95% bootstrap intervals (1000 replicates over subjects). Cox on Tier B and elastic net stay near 1; the untuned random survival forest over-predicts. Arrows mark intervals cut at the axis.\n")
  cat("Fig 2. O/E by protocol era for Cox-B without and with an era term, 2-year horizon. ADNI-GO (3 events by 2 years at baseline) and ADNI4 (2-year status observable for few subjects; estimate unstable) are omitted. Numbers at the bottom are observed events. Models with an era term were fitted at baseline and 12 months only. Points based on fewer than 5 events are not drawn.\n")
  cat("Fig 3. Subgroup O/E at baseline, 2 years, over the grey range that ADNI could not distinguish from 1 with 80% power after allowing for subjects lost before 2 years. Estimates for Black and Hispanic participants lie entirely within it.\n")
  cat("Fig 4. Expected events among subjects with known 2-year status (dots) and observed events (crosses) per subgroup at baseline, against the expected events needed to detect a given O/E with 80% power. The power calculation treats expected events as fixed and uses a Poisson approximation.\n")
  cat("Fig 5. Net benefit above the better of treat-all and treat-none at 2 years, by landmark and risk threshold. The shaded band is the 95% interval for Cox-B.\n")
  cat("Fig 6. Coverage of the true outcome class by class-conditional cross-conformal sets (target 90%), censoring-weighted, by subgroup at baseline for dementia within 2 years. Hollow points have fewer than 30 events. Intervals reflect evaluation resampling only.\n")
  cat("Fig 7. Change in out-of-fold AUC relative to Cox-B fitted on the same rows with plasma available (209 baseline rows, 62 at 12 months), with paired bootstrap intervals. The shuffled-plasma row is a negative control. Panels with fewer than 10 events (12 months, 1 year: 2 events) are not shown.\n")
  cat("Fig S1. As Fig 6 for marginal conformal sets, which reach 90% overall but cover only about 54% of events.\n")
  cat("Fig S2. Difference between 1 - Kaplan-Meier risk (death censored) and the Aalen-Johansen cumulative incidence (death competing), whole cohort.\n")
})

section("4. Checks", {
  ok_in <- !is.null(oe) && !is.null(dc) && !is.null(cf) && !is.null(pg) && !is.null(cr) && !is.null(p11)
  check("all input CSVs loaded", ok_in)
  if (ok_in) {
    a <- oe[oe$variable == "all" & oe$model == "Cox-B" & oe$landmark == "bl" & oe$horizon == 2, ]
    check("whole-cohort baseline events at 2 years are 201", nrow(a) == 1 && a$events_by_h == 201, sprintf("(%s)", paste(a$events_by_h, collapse = ",")))
    e3 <- oe[oe$variable == "origprot" & oe$level == "ADNI3" & oe$model == "Cox-B" & oe$landmark == "bl" & oe$horizon == 2, ]
    check("ADNI3 O/E for Cox-B at baseline, 2 years matches the 08 report (0.502)", nrow(e3) == 1 && abs(e3$O_over_E - 0.502) < 0.005, sprintf("(%s)", paste(e3$O_over_E, collapse = ",")))
    nb <- dc[dc$model == "Cox-B" & dc$landmark == "bl" & dc$horizon == 2 & abs(dc$threshold - 0.20) < 1e-9, ]
    check("Cox-B net benefit at baseline, 2 years, threshold 0.20 matches the 09 report (0.103)", nrow(nb) == 1 && abs(nb$NB_model - 0.103) < 0.002, sprintf("(%s)", paste(nb$NB_model, collapse = ",")))
    bb <- cf[cf$landmark == "bl" & cf$model == "Cox-B" & cf$horizon == 2 & cf$method == "class-cond-IPCW" & cf$variable == "race_bw" & cf$level == "Black", ]
    check("class-conditional coverage for Black participants matches the 09 report (0.83)", nrow(bb) == 1 && abs(bb$coverage_ipcw - 0.83) < 0.005, sprintf("(%s)", paste(bb$coverage_ipcw, collapse = ",")))
    pl <- pg[pg$landmark == "bl" & pg$model == "Cox-B+plasma" & pg$horizon == 2, ]
    check("plasma AUC gain at baseline, 2 years matches the 10 report (0.019)", nrow(pl) == 1 && abs(pl$delta_AUC_vs_CoxB - 0.019) < 0.002, sprintf("(%s)", paste(pl$delta_AUC_vs_CoxB, collapse = ",")))
    check("whole-cohort KM minus Aalen-Johansen never exceeds 0.4 points", max(cr$difference_pp[cr$variable == "all"], na.rm = TRUE) < 0.4,
          sprintf("(max %.3f)", max(cr$difference_pp[cr$variable == "all"], na.rm = TRUE)))
    m2 <- max(abs(cr$difference_pp[cr$horizon == 2]), na.rm = TRUE)
    check("no subgroup differs by 1 point at 2 years (10 report)", m2 < 1, sprintf("(max %.3f)", m2))
    check("every model drawn in Fig 1 is present", all(M4 %in% oe$model))
    check("every era drawn in Fig 2 is present at baseline", all(ERAS %in% oe$level[oe$variable == "origprot" & oe$landmark == "bl"]))
  }
  fdir <- file.path(CFG$dir_out, "figures")
  for (nm in names(fig_ok)) {
    p <- file.path(fdir, paste0(nm, c(".png", ".pdf")))
    check(paste("figure files written and non-trivial:", nm), isTRUE(fig_ok[[nm]]) && all(file.exists(p)) && all(file.size(p) > 5000))
  }
  cat(sprintf("\n%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  utils::write.csv(tests, A("12_tests.csv"), row.names = FALSE)
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
