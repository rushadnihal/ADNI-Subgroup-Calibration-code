# 03_outcomes.R ------------------------------------------------------------
# Builds the landmark cohort and outcomes from output/interim/visits.rds.
# One row per subject x landmark (bl, m12, m24, m36) at which the subject is
# diagnosed MCI. Time zero is the landmark visit date.
#
# Outcome: first dementia diagnosis after the landmark date.
#   status 1 = dementia, 2 = death (competing), 0 = censored at last diagnosis visit.
# Death counts as a competing event only if it falls within DEATH_WINDOW_DAYS of
# the last diagnosis visit; later deaths are censored at the last visit, because
# dementia status between the last visit and a late death is unknown.
# A second version (status_cens) treats death as plain censoring (comparator).
# A third (status_conf) uses confirmed dementia (next visit also dementia, or last visit).
#
# Writes: output/interim/outcomes.rds   (has subject IDs, stays local)
#         output/audit/03_*.csv, 03_outcomes_report.txt   (counts only)
# Run:    source("R/03_outcomes.R")   (after 02_build_visits.R; needs survival, which ships with R)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")
if (!requireNamespace("survival", quietly = TRUE))
  stop("Package 'survival' is missing. It ships with R; run install.packages('survival').")

report <- file.path(CFG$dir_out, "audit", "03_outcomes_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Outcomes build | run:", format(Sys.time()), "\n")

# ---- settings you may want to change --------------------------------------
DAYS              <- 365.25
GRACE_YEARS       <- 0.25   # a subject followed to within 3 months of the horizon counts as followed to it
DEATH_WINDOW_DAYS <- 365    # death is a competing event only within this long of the last diagnosis visit
HORIZONS          <- CFG$horizons_years
SUB_VARS          <- c("sex", "race", "hisp", "edu_grp", "apoe4_grp", "age_grp", "origprot")
# ---------------------------------------------------------------------------

visits <- readRDS(file.path(CFG$dir_out, "interim", "visits.rds"))
need <- c("RID", "EXAMDATE", "DIAGNOSIS", "VISCODE2", "origprot", "sex", "race", "hisp",
          "edu_grp", "apoe4_grp", "age_visit", "died", "death_date")
stopifnot(all(need %in% names(visits)))
data_cutoff <- max(visits$EXAMDATE, na.rm = TRUE)
cat("Data cut-off (latest diagnosis visit):", format(data_cutoff), "\n")
cat("ORIGPROT values:", paste(names(table(visits$origprot, useNA = "ifany")), collapse = ", "), "\n")

oc <- NULL
section("1. Landmark cohort and flow", {
  subj <- visits %>% arrange(RID, EXAMDATE) %>% group_by(RID) %>%
    summarise(last_obs = max(EXAMDATE), n_visits = n(),
              died = first(died), death_date = first(death_date)) %>% ungroup() %>%
    mutate(death_known = died & !is.na(death_date),
           death_inconsistent = death_known & death_date < last_obs,
           death_eff = as.Date(ifelse(!death_known, NA_real_,
                                      pmax(as.numeric(death_date), as.numeric(last_obs))),
                               origin = "1970-01-01"))

  lm <- visits %>% filter(VISCODE2 %in% CFG$landmarks, DIAGNOSIS == 2) %>%
    group_by(RID, VISCODE2) %>% slice(1) %>% ungroup() %>%
    transmute(RID, landmark = factor(VISCODE2, levels = CFG$landmarks), t0 = EXAMDATE,
              origprot, sex, race, hisp, edu_grp, apoe4_grp, age_lm = age_visit)
  lm$age_grp <- as.character(cut(lm$age_lm, c(-Inf, 65, 75, Inf), right = FALSE,
                                 labels = c("<65", "65-74", "75+")))

  j <- lm %>% select(RID, landmark, t0) %>%
    inner_join(visits %>% select(RID, EXAMDATE, DIAGNOSIS), by = "RID", relationship = "many-to-many")
  prior <- j %>% filter(EXAMDATE <= t0, DIAGNOSIS == 3) %>% distinct(RID, landmark) %>%
    mutate(prior_dementia = TRUE)
  after <- j %>% filter(EXAMDATE > t0) %>% arrange(RID, landmark, EXAMDATE) %>%
    group_by(RID, landmark) %>% mutate(next_dx = lead(DIAGNOSIS)) %>% ungroup()
  first_dem <- after %>% filter(DIAGNOSIS == 3) %>% group_by(RID, landmark) %>%
    summarise(t_dem = min(EXAMDATE)) %>% ungroup()
  conf_dem <- after %>% filter(DIAGNOSIS == 3, is.na(next_dx) | next_dx == 3) %>%
    group_by(RID, landmark) %>% summarise(t_dem_conf = min(EXAMDATE)) %>% ungroup()

  oc <- lm %>% left_join(subj, by = "RID") %>%
    left_join(first_dem, by = c("RID", "landmark")) %>%
    left_join(conf_dem,  by = c("RID", "landmark")) %>%
    left_join(prior,     by = c("RID", "landmark")) %>%
    mutate(prior_dementia = coalesce(prior_dementia, FALSE),
           has_followup = last_obs > t0,
           in_set = !prior_dementia & has_followup,   # modelling set
           in_audit = !prior_dementia)                # censoring audit keeps zero-follow-up rows

  build <- function(tdem) {
    dem <- !is.na(tdem)
    gap <- as.numeric(oc$death_eff - oc$last_obs)
    dth <- !dem & !is.na(gap) & gap <= DEATH_WINDOW_DAYS
    end <- ifelse(dem, as.numeric(tdem),
                  ifelse(dth, as.numeric(oc$death_eff), as.numeric(oc$last_obs)))
    list(time = (end - as.numeric(oc$t0)) / DAYS,
         status = ifelse(dem, 1L, ifelse(dth, 2L, 0L)))
  }
  p <- build(oc$t_dem);      oc$time <- p$time;           oc$status <- p$status
  p <- build(oc$t_dem_conf); oc$time_conf <- p$time;      oc$status_conf <- p$status
  oc$time_cens   <- (ifelse(!is.na(oc$t_dem), as.numeric(oc$t_dem), as.numeric(oc$last_obs)) -
                       as.numeric(oc$t0)) / DAYS
  oc$status_cens <- as.integer(!is.na(oc$t_dem))

  for (h in HORIZONS) {
    obs <- (as.numeric(oc$t0) + h * DAYS) <= as.numeric(data_cutoff)
    oc[[paste0("obs_h", h)]] <- obs
    oc[[paste0("class_h", h)]] <-
      ifelse(oc$status == 1 & oc$time <= h, "dementia",
      ifelse(oc$status == 2 & oc$time <= h, "death",
      ifelse(oc$time >= h - GRACE_YEARS, "event_free",
      ifelse(!obs, "administrative", "lost"))))
  }

  flow <- oc %>% group_by(landmark) %>% summarise(
    n_mci_at_landmark   = n(),
    n_prior_dementia    = sum(prior_dementia),
    n_no_followup       = sum(!prior_dementia & !has_followup),
    n_analysis_set      = sum(in_set),
    n_first_dem_unconfirmed = sum(in_set & !is.na(t_dem) &
                                    (is.na(t_dem_conf) | t_dem_conf != t_dem), na.rm = TRUE),
    n_death_competing   = sum(in_set & status == 2),
    n_death_date_before_last_visit = sum(in_set & death_inconsistent)) %>% ungroup()
  show(flow); save_csv(flow, "03_cohort_flow.csv")
  cat("\nAnalysis set: subject-landmark rows =", sum(oc$in_set),
      "| unique subjects =", n_distinct(oc$RID[oc$in_set]), "\n")
})

section("2. Events and follow-up by landmark and horizon (audit set: includes subjects with no visit after the landmark)", {
  d <- oc[oc$in_audit, ]
  ev <- do.call(rbind, lapply(HORIZONS, function(h) {
    do.call(rbind, lapply(CFG$landmarks, function(l) {
      x <- d[d$landmark == l, ]; cl <- x[[paste0("class_h", h)]]
      data.frame(landmark = l, horizon = h, n = nrow(x),
                 dementia = sum(cl == "dementia"),
                 dementia_confirmed = sum(x$status_conf == 1 & x$time_conf <= h),
                 death_competing = sum(cl == "death"),
                 event_free_at_h = sum(cl == "event_free"),
                 lost_before_h = sum(cl == "lost"),
                 of_which_no_visit_after_landmark = sum(cl == "lost" & !x$has_followup),
                 administrative = sum(cl == "administrative"))
    }))
  }))
  ev$pct_lost_of_observable <- round(100 * ev$lost_before_h / (ev$n - ev$administrative), 1)
  ev$pct_lost_excl_no_visit <- round(100 * (ev$lost_before_h - ev$of_which_no_visit_after_landmark) /
                                       (ev$n - ev$administrative - ev$of_which_no_visit_after_landmark), 1)
  show(ev, 40); save_csv(ev, "03_events_by_horizon.csv")
  cat("\nReference (01_audit, before 6 conflicting baseline rows were dropped): bl dementia 80/202/278; m36 16/41/51.\n")
  cat("Classes: dementia/death = event by h; event_free_at_h = followed to within", GRACE_YEARS,
      "y of h without event; lost_before_h = observable but follow-up ended earlier;\n",
      "administrative = h not yet reachable before the data cut-off.\n")
})

section("3. Death dates and the competing-event window", {
  s <- oc %>% filter(in_set) %>% distinct(RID, .keep_all = TRUE)
  gap <- as.numeric(s$death_date - s$last_obs)
  b <- cut(gap, c(-Inf, -1, 180, 365, 730, Inf),
           labels = c("death date BEFORE last diagnosis visit", "0-180 d after", "181-365 d after",
                      "366-730 d after", ">730 d after"))
  tb <- data.frame(group = c("flagged died, no death date", levels(b), "not recorded as dead"),
                   n = c(sum(s$died & is.na(s$death_date)), as.integer(table(b)), sum(!s$died)))
  tb$used_as_competing_event <- c("no (censored)", "yes (time moved to last visit)", "yes", "yes",
                                  "no (censored at last visit)", "no (censored at last visit)", "n/a")
  show(tb); save_csv(tb, "03_death_gap.csv")
  cat("\nsubjects in analysis set:", nrow(s), "(any landmark). Change DEATH_WINDOW_DAYS at the top to move the cut.\n")
})

long <- NULL
section("4. Who is lost to follow-up before the horizon? (observable subjects only)", {
  d <- oc[oc$in_audit, ]
  long <- do.call(rbind, lapply(HORIZONS, function(h) {
    x <- d[, c("RID", "landmark", "has_followup", SUB_VARS)]
    x$horizon <- h; x$class <- d[[paste0("class_h", h)]]; x
  }))
  sub_tab <- do.call(rbind, lapply(SUB_VARS, function(v) {
    long %>% filter(class != "administrative") %>%
      group_by(landmark, horizon, level = as.character(.data[[v]])) %>%
      summarise(n_observable = n(), dementia = sum(class == "dementia"),
                death = sum(class == "death"), event_free = sum(class == "event_free"),
                lost = sum(class == "lost"),
                lost_no_visit_after_landmark = sum(class == "lost" & !has_followup)) %>% ungroup() %>%
      mutate(variable = v, pct_lost = round(100 * lost / n_observable, 1),
             pct_lost_excl_no_visit = round(100 * (lost - lost_no_visit_after_landmark) /
                                              (n_observable - lost_no_visit_after_landmark), 1))
  }))
  show(sub_tab %>% filter(landmark == "bl", horizon == 2), 60)
  save_csv(sub_tab, "03_lost_by_subgroup.csv")

  rp <- long %>% filter(class != "administrative") %>%
    mutate(race2 = ifelse(race == "Black", "Black", "Other")) %>%
    group_by(landmark, horizon, origprot, race2) %>%
    summarise(n_observable = n(), lost = sum(class == "lost"),
              lost_no_visit_after_landmark = sum(class == "lost" & !has_followup)) %>% ungroup() %>%
    mutate(pct_lost = round(100 * lost / n_observable, 1))
  cat("\nBlack vs other by ORIGPROT (bl landmark, 2-year horizon):\n")
  show(rp %>% filter(landmark == "bl", horizon == 2), 40)
  save_csv(rp, "03_lost_black_vs_other_by_origprot.csv")
})

section("5. Reverse Kaplan-Meier: probability of being lost before the horizon (observable, ADNI4 excluded)", {
  rkm <- function(l, h, v) {
    d <- oc[oc$in_audit & oc$landmark == l & !(oc$origprot %in% "ADNI4") & oc[[paste0("obs_h", h)]], ]
    d$lvl <- if (v == "race_bw") ifelse(d$race %in% c("Black", "White"), d$race, NA_character_)
             else as.character(d[[v]])
    d <- d[!is.na(d$lvl) & !(d$lvl %in% c("Unknown", "missing")), ]   # unknown groups would distort the log-rank test
    d$tt <- pmin(d$time, h)
    d$ev <- as.integer(d$status == 0 & d$time < h - GRACE_YEARS)
    if (length(unique(d$lvl)) < 2) return(NULL)
    f  <- survival::survfit(survival::Surv(tt, ev) ~ lvl, data = d)
    s  <- summary(f, times = h, extend = TRUE)
    sd <- survival::survdiff(survival::Surv(tt, ev) ~ lvl, data = d)
    lv <- sub("^lvl=", "", as.character(s$strata))
    data.frame(landmark = l, horizon = h, variable = v, level = lv,
               n = as.integer(table(d$lvl)[lv]),
               pct_lost_km = round(100 * (1 - s$surv), 1),
               lcl = round(100 * (1 - s$upper), 1), ucl = round(100 * (1 - s$lower), 1),
               logrank_p = signif(1 - pchisq(sd$chisq, length(sd$n) - 1), 3))
  }
  grid <- expand.grid(l = c("bl", "m12"), h = HORIZONS[HORIZONS >= 2],
                      v = c(setdiff(SUB_VARS, "origprot"), "race_bw"), stringsAsFactors = FALSE)
  out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
    tryCatch(rkm(grid$l[i], grid$h[i], grid$v[i]), error = function(e) NULL)))
  show(out %>% filter(landmark == "bl", horizon == 2), 40)
  save_csv(out, "03_reverse_km_by_subgroup.csv")
  cat("\nSmall cells give wide intervals; read the p-values with the counts beside them.\n")
})

section("6. Save", {
  saveRDS(oc, file.path(CFG$dir_out, "interim", "outcomes.rds"))
  cat("Saved output/interim/outcomes.rds:", nrow(oc), "rows,", ncol(oc), "columns (in_set marks the analysis set).\n")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
