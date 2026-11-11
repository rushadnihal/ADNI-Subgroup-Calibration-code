# 11_events_needed.R ---------------------------------------------------------
# How many events does a subgroup calibration audit need, and what does ADNI have?
#
# Uses the Cox-B out-of-fold predictions saved by 07 (no models are refit).
#
# 1. Events needed. For calibration-in-the-large (observed/expected, O/E), a subgroup
#    with E expected events can detect an O/E of m with 80% power (two-sided 5% test,
#    Poisson approximation) when
#        E >= (z_{0.975} + z_{0.80} * sqrt(m))^2 / (m - 1)^2.
#    Tabulated for a grid of m, and converted to participants needed at several risks.
# 2. ADNI against that requirement, for every landmark and subgroup level at the 2-year
#    horizon: subjects, observed events, expected events, the share of subjects whose
#    status by the horizon is known (not lost earlier), and the O/E range each group could
#    detect with 80% power, with and without allowing for loss to follow-up.
# 3. The same for contrasts between two groups (ratio of O/E, as in 08): the smallest
#    ratio detectable with 80% power, set against the intervals 08 found.
#
# The power figures ignore model error and treat expected events as fixed, so they are
# planning guides, not exact. Death is treated as censoring, as in 04 to 10.
#
# Writes: output/audit/11_*.csv, 11_events_needed_report.txt
#         output/figures/11_events_vs_needed.png (counts only)
# Run:    source("R/11_events_needed.R")   (after 07 and 08; under a minute)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
K <- 5; HORIZONS <- CFG$horizons_years
source(if (file.exists("R/utils_eval.R")) "R/utils_eval.R" else "utils_eval.R")

report <- file.path(CFG$dir_out, "audit", "11_events_needed_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Events needed for subgroup calibration audits | run:", format(Sys.time()), "\n")

# ---- settings -------------------------------------------------------------
MODEL   <- "Cox-B"
H       <- 2                       # horizon (years) for the subgroup table
MIN_N   <- 10                      # list levels down to this many subjects
ZA      <- qnorm(0.975); ZB <- qnorm(0.80)
M_GRID  <- c(0.5, 0.6, 0.67, 0.75, 0.8, 0.9, 1.1, 1.25, 1.33, 1.5, 2)
RISKS   <- c(0.10, 0.20, 0.30)     # 2-year risks used to convert events to participants
VARS    <- c("sex", "race_bw", "hisp", "edu_grp", "apoe4_grp", "age_grp", "origprot")
CONTRASTS <- list(c("race_bw", "Black", "White"), c("hisp", "Hispanic", "Not Hispanic"),
                  c("sex", "Female", "Male"), c("edu_grp", "<=12y", ">12y"),
                  c("apoe4_grp", "1+", "0"), c("age_grp", "75+", "<65"),
                  c("origprot", "ADNI3", "ADNI1"), c("origprot", "ADNI4", "ADNI1"))
# ---------------------------------------------------------------------------

tests <- data.frame(test = character(), result = character(), detail = character(), stringsAsFactors = FALSE)
check <- function(name, ok, detail = "") {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s %s\n", if (ok) "PASS" else "FAIL", name, detail))
  tests <<- rbind(tests, data.frame(test = name, result = if (ok) "PASS" else "FAIL", detail = detail, stringsAsFactors = FALSE))
  invisible(ok)
}
rd <- function(x) { nn <- vapply(x, is.numeric, logical(1)); x[nn] <- lapply(x[nn], round, 3); x }

# expected events needed to detect O/E = m with 80% power
E_needed <- function(m) (ZA + ZB * sqrt(m))^2 / (m - 1)^2
# closest O/E to 1 on each side detectable with expected events E (NA if none in range)
detect_range <- function(E) {
  if (!is.finite(E) || E <= 0) return(c(below = NA_real_, above = NA_real_))
  lo <- seq(0.995, 0.01, by = -0.005); hi <- seq(1.005, 10, by = 0.005)
  b <- lo[which(E_needed(lo) <= E)[1]]; a <- hi[which(E_needed(hi) <= E)[1]]
  c(below = if (length(b)) b else NA_real_, above = if (length(a)) a else NA_real_)
}
# smallest ratio of two O/E values detectable with 80% power, expected events E1 and E2
detect_ratio <- function(E1, E2) {
  if (!all(is.finite(c(E1, E2))) || E1 <= 0 || E2 <= 0) return(c(below = NA_real_, above = NA_real_))
  d <- (ZA + ZB) * sqrt(1 / E1 + 1 / E2); c(below = exp(-d), above = exp(d))
}

o07 <- readRDS(file.path(CFG$dir_out, "interim", "oof_predictions.rds"))
D <- o07$d; Pm <- o07$P; stopifnot(identical(as.numeric(o07$horizons), as.numeric(HORIZONS)), MODEL %in% names(Pm[["bl"]]))
LM <- CFG$landmarks; j <- match(H, HORIZONS)

section("1. Events needed to detect a given O/E (80% power, two-sided 5% test)", {
  grid <- data.frame(O_over_E_to_detect = M_GRID, expected_events_needed = ceiling(E_needed(M_GRID)))
  for (r in RISKS) grid[[paste0("participants_needed_at_risk_", r)]] <- ceiling(grid$expected_events_needed / r)
  save_csv(grid, "11_events_needed_grid.csv"); assign("grid", grid, envir = .GlobalEnv)
  show(grid, 20)
  cat("\nReading: a 20% miscalibration (O/E of 0.8 or 1.25) needs about",
      ceiling(E_needed(0.8)), "or", ceiling(E_needed(1.25)),
      "expected events in the subgroup. At a 20% two-year risk that is roughly",
      ceiling(E_needed(0.8) / 0.2), "or", ceiling(E_needed(1.25) / 0.2), "participants with complete follow-up per group.\n")
})

sub_rows <- list(); con_rows <- list()
section("2. What ADNI has: subgroup levels at the 2-year horizon", {
  for (l in LM) {
    dl <- D[D$landmark == l, ]; n <- nrow(dl); p <- Pm[[l]][[MODEL]][, j]
    known <- (dl$ev == 1 & dl$time <= H) | dl$time > H
    levels_l <- list(list(var = "all", level = "all", idx = seq_len(n)))
    for (v in VARS) {
      lv <- as.character(dl[[v]])
      for (u in setdiff(unique(lv), c(NA, "Unknown", "missing"))) {
        i <- which(lv == u); if (length(i) >= MIN_N) levels_l[[length(levels_l) + 1]] <- list(var = v, level = u, idx = i)
      }
    }
    for (s in levels_l) {
      i <- s$idx; ev <- sum(dl$ev[i] == 1 & dl$time[i] <= H)
      E_all <- sum(p[i]); E_eff <- sum(p[i][known[i]])
      d1 <- detect_range(E_all); d2 <- detect_range(E_eff)
      sub_rows[[length(sub_rows) + 1]] <- data.frame(landmark = l, variable = s$var, level = s$level, n = length(i),
        events_by_h = ev, share_status_known = mean(known[i]), expected_events_all = E_all, expected_events_known = E_eff,
        needed_for_OE_0.8 = E_needed(0.8), needed_for_OE_1.25 = E_needed(1.25),
        fraction_of_needed_0.8 = E_eff / E_needed(0.8),
        detectable_OE_below_ideal = d1[["below"]], detectable_OE_above_ideal = d1[["above"]],
        detectable_OE_below_after_loss = d2[["below"]], detectable_OE_above_after_loss = d2[["above"]],
        stringsAsFactors = FALSE)
    }
    for (ct in CONTRASTS) {
      a <- Filter(function(s) s$var == ct[1] && s$level == ct[2], levels_l); b <- Filter(function(s) s$var == ct[1] && s$level == ct[3], levels_l)
      if (length(a) != 1 || length(b) != 1) next
      ia <- a[[1]]$idx; ib <- b[[1]]$idx
      E1 <- sum(p[ia][known[ia]]); E2 <- sum(p[ib][known[ib]])
      dr <- detect_ratio(E1, E2)
      con_rows[[length(con_rows) + 1]] <- data.frame(landmark = l, variable = ct[1], group = ct[2], reference = ct[3],
        expected_events_group = E1, expected_events_reference = E2,
        detectable_ratio_below = dr[["below"]], detectable_ratio_above = dr[["above"]], stringsAsFactors = FALSE)
    }
  }
  sub_df <- rd(do.call(rbind, sub_rows)); con_df <- rd(do.call(rbind, con_rows))
  save_csv(sub_df, "11_subgroup_events_vs_needed.csv"); save_csv(con_df, "11_contrast_detectable_ratio.csv")
  assign("sub_df", sub_df, envir = .GlobalEnv); assign("con_df", con_df, envir = .GlobalEnv)
  cat("Baseline, 2 years (events_by_h observed; expected_events_known = predicted events among subjects whose status is known):\n")
  show(sub_df[sub_df$landmark == "bl", c("variable", "level", "n", "events_by_h", "share_status_known", "expected_events_known",
        "fraction_of_needed_0.8", "detectable_OE_below_after_loss", "detectable_OE_above_after_loss")], 60)
  ok <- sub_df$expected_events_known >= E_needed(0.8) & sub_df$expected_events_known >= E_needed(1.25)
  cat(sprintf("\nSubgroup levels (all landmarks, 2 years, listed down to %d subjects) with enough events to detect a 20%% miscalibration either way: %d of %d.\n",
              MIN_N, sum(ok), nrow(sub_df)))
  cat("The levels that do:\n"); show(sub_df[ok, c("landmark", "variable", "level", "n", "events_by_h", "expected_events_known")], 40)
})

section("3. Contrasts between groups: smallest detectable ratio of O/E", {
  cat("Baseline, 2 years. Compare with the intervals in 08_oe_contrasts.csv:\n")
  show(con_df[con_df$landmark == "bl", ], 20)
  r8 <- tryCatch(readr::read_csv(file.path(CFG$dir_out, "audit", "08_oe_contrasts.csv"), show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(r8)) {
    x <- r8[r8$model == MODEL & r8$horizon == H & r8$landmark == "bl" & r8$variable %in% c("race_bw", "hisp"), c("variable", "group", "reference", "ratio_of_OE", "lcl", "ucl")]
    cat("\nObserved ratios of O/E at baseline (from 08):\n"); show(x, 10)
  }
})

section("4. Figure", {
  fdir <- file.path(CFG$dir_out, "figures"); dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
  ok1 <- tryCatch({
    d <- sub_df[sub_df$landmark == "bl", ]; d <- d[order(d$expected_events_known), ]
    lab <- ifelse(d$variable == "all", "all", paste0(d$variable, ": ", d$level))
    grDevices::png(file.path(fdir, "11_events_vs_needed.png"), width = 1500, height = 1100, res = 150)
    par(mar = c(4, 11, 3, 1))
    yy <- seq_len(nrow(d))
    plot(NA, xlim = c(1, max(d$expected_events_known, E_needed(0.8)) * 1.3), ylim = c(0.5, nrow(d) + 0.5), log = "x", yaxt = "n",
         xlab = "Expected events among subjects with known status (log scale)", ylab = "",
         main = "Baseline, 2 years: events available vs events needed to detect miscalibration", cex.main = 0.85)
    axis(2, at = yy, labels = lab, las = 1, cex.axis = 0.65)
    segments(1, yy, d$expected_events_known, yy, col = "#3b6ea5", lwd = 5, lend = 1)
    points(d$events_by_h, yy, pch = 4, col = "#c0504d", lwd = 1.5)
    for (m in c(0.67, 0.8)) abline(v = E_needed(m), lty = 2, col = "grey40")
    mtext(sprintf("Dashed lines: events needed for O/E = 0.67 (%d) and 0.8 (%d)", ceiling(E_needed(0.67)), ceiling(E_needed(0.8))), side = 3, line = 0, cex = 0.7)
    legend("bottomright", legend = c("Expected events, known status", "Observed events by 2 years"), pch = c(15, 4), col = c("#3b6ea5", "#c0504d"), bty = "n", cex = 0.7)
    grDevices::dev.off(); TRUE
  }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); cat("  figure not written:", conditionMessage(e), "\n"); FALSE })
  cat("Figure written to output/figures/ :", if (ok1) "11_events_vs_needed.png" else "none", "\n")
})

section("5. Checks", {
  set.seed(CFG$seed + 1300)
  Eneed <- E_needed(0.8); o <- rpois(20000, 0.8 * Eneed)
  pw <- mean(abs(o - Eneed) / sqrt(Eneed) > ZA)
  check("simulation: the events formula gives about 80% power (0.75 to 0.85)", pw > 0.75 && pw < 0.85, sprintf("(power %.3f at E = %.0f)", pw, Eneed))
  o0 <- rpois(20000, 200); t1 <- mean(abs(o0 - 200) / sqrt(200) > ZA)
  check("simulation: the test rejects about 5% of the time when calibrated (0.03 to 0.07)", t1 > 0.03 && t1 < 0.07, sprintf("(rate %.3f)", t1))
  E1 <- 100; E2 <- 100; r <- detect_ratio(E1, E2)[["above"]]
  o1 <- rpois(20000, r * E1); o2 <- rpois(20000, E2); ok <- o1 > 0 & o2 > 0
  pw2 <- mean(abs(log((o1[ok] / E1) / (o2[ok] / E2))) / sqrt(1 / o1[ok] + 1 / o2[ok]) > ZA)
  check("simulation: the contrast formula gives about 80% power (0.72 to 0.88)", pw2 > 0.72 && pw2 < 0.88, sprintf("(power %.3f at ratio %.2f)", pw2, r))
  Es <- c(20, 50, 100, 200, 400); dr <- sapply(Es, function(e) detect_range(e)[["above"]])
  check("the detectable range narrows as expected events grow", all(diff(dr) < 0))
  r8 <- tryCatch(readr::read_csv(file.path(CFG$dir_out, "audit", "08_oe_over_time.csv"), show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(r8)) {
    a <- r8[r8$variable == "all" & r8$model == MODEL & r8$landmark == "bl" & r8$horizon == H, ]
    s <- sub_df[sub_df$landmark == "bl" & sub_df$variable == "all", ]
    check("whole-cohort baseline events at 2 years match 08", nrow(a) == 1 && nrow(s) == 1 && a$events_by_h == s$events_by_h, sprintf("(%s vs %s)", a$events_by_h, s$events_by_h))
    dl <- D[D$landmark == "bl", ]; mp <- mean(Pm[["bl"]][[MODEL]][, j])
    check("mean predicted risk reproduces 08 (within 0.002)", abs(mp - a$mean_pred) < 0.002, sprintf("(%.4f vs %.4f)", mp, a$mean_pred))
  } else check("08 results found for the regression check", FALSE)
  check("saved predictions line up with the outcome rows at every landmark",
        all(sapply(LM, function(l) nrow(D[D$landmark == l, ]) == nrow(Pm[[l]][[1]]))))
  cat(sprintf("\n%d of %d checks passed.\n", sum(tests$result == "PASS"), nrow(tests)))
  save_csv(tests, "11_tests.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
