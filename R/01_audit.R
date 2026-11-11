# 01_audit.R ---------------------------------------------------------------
# Read-only data audit. No modelling. Writes aggregate tables (no subject IDs)
# to output/audit/ and a full text log to output/audit/01_audit_report.txt.
#
# Run:  setwd to analysis1111, then  source("R/01_audit.R")
# Takes about 1-3 minutes (UCSFFSX7 is 33 MB).
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")

report <- file.path(CFG$dir_out, "audit", "01_audit_report.txt")
sink(report, split = TRUE)
cat("ADNI data audit | run:", format(Sys.time()), "\n")
cat(R.version.string, "\n")
cat("dplyr", as.character(packageVersion("dplyr")),
    "| readr", as.character(packageVersion("readr")), "\n")
options(width = 170)

# ---- A. file inventory ----------------------------------------------------
section("A. File inventory", {
  f   <- list.files(CFG$dir_data, full.names = TRUE)
  inv <- data.frame(file = basename(f), mb = round(file.size(f) / 1e6, 3),
                    is_dir = dir.exists(f), stringsAsFactors = FALSE)
  show(inv); save_csv(inv, "01_file_inventory.csv")
  rda  <- sub("\\.rda$", "", list.files(CFG$dir_rda, pattern = "\\.rda$"))
  need <- c("DXSUM", "PTDEMOG", "APOERES", "DM", "REGISTRY", "CDR", "ADAS", "FAQ", "MMSE")
  cat("\nADNIMERGE2 .rda tables found:", length(rda), "\n")
  cat("Needed tables present:\n")
  print(data.frame(table = need, present = need %in% rda), row.names = FALSE)
})

# ---- B. load tables -------------------------------------------------------
tb <- list()
load_try <- function(name, expr) {
  r <- tryCatch(expr, error = function(e) {
    cat("  FAILED ", name, ": ", conditionMessage(e), "\n", sep = ""); NULL })
  if (!is.null(r)) {
    tb[[name]] <<- r
    cat(sprintf("  %-10s %8d rows x %4d cols\n", name, nrow(r), ncol(r)))
  }
}
section("B. Loading tables", {
  load_try("DXSUM",    read_adni_csv("DXSUM_11Nov2026.csv"))
  load_try("PLASMA",   read_adni_csv("UPENN_PLASMA_FUJIREBIO_QUANTERIX_11Nov2026.csv"))
  load_try("AV45",     read_adni_csv("BAIPETNMRCAV45_10_23_20_11Nov2026.csv"))
  load_try("NEUROBAT", read_adni_csv("NEUROBAT_11Nov2026.csv"))
  load_try("UCSFFSX7", read_adni_csv("UCSFFSX7_11Nov2026.csv"))
  load_try("DATADIC",  read_adni_csv("DATADIC_11Nov2026.csv"))
  for (nm in c("PTDEMOG", "APOERES", "DM", "REGISTRY", "CDR", "ADAS", "FAQ", "MMSE")) {
    load_try(nm, load_rda(nm))
  }
})

# ---- C. per-table summary -------------------------------------------------
section("C. Table summaries (keys, duplicates, dates, phases)", {
  specs <- list(
    list("DXSUM",    "RID",    "VISCODE2", "EXAMDATE", "PHASE"),
    list("PLASMA",   "RID",    "VISCODE2", "EXAMDATE", "PHASE"),
    list("AV45",     "RID",    "VISCODE2", "EXAMDATE", "COLPROT"),
    list("NEUROBAT", "RID",    "VISCODE2", "VISDATE",  "PHASE"),
    list("UCSFFSX7", "RID",    "VISCODE2", "EXAMDATE", "PHASE"),
    list("PTDEMOG",  "RID",    "VISCODE2", "VISDATE",  "COLPROT"),
    list("APOERES",  "RID",    "VISCODE",  "APTESTDT", "COLPROT"),
    list("REGISTRY", "RID",    "VISCODE2", "EXAMDATE", "COLPROT"),
    list("CDR",      "RID",    "VISCODE2", "VISDATE",  "COLPROT"),
    list("ADAS",     "RID",    "VISCODE2", "VISDATE",  "COLPROT"),
    list("FAQ",      "RID",    "VISCODE2", "VISDATE",  "COLPROT"),
    list("MMSE",     "RID",    "VISCODE2", "VISDATE",  "COLPROT"),
    list("DM",       "SUBJID", NA,         "DMDTC",    "ORIGPROT")
  )
  summ_one <- function(s) {
    df <- tb[[s[[1]]]]
    if (is.null(df)) return(NULL)
    id <- s[[2]]; vis <- s[[3]]; dt <- s[[4]]; ph <- s[[5]]
    has <- function(cn) !is.na(cn) && cn %in% names(df)
    d  <- if (has(dt)) to_date(df[[dt]]) else NULL
    tt <- if (has(ph)) table(df[[ph]], useNA = "ifany") else NULL
    data.frame(
      table        = s[[1]],
      rows         = nrow(df),
      cols         = ncol(df),
      subjects     = if (has(id)) length(unique(df[[id]])) else NA,
      dup_id_visit = if (has(id) && has(vis)) sum(duplicated(df[c(id, vis)])) else NA,
      na_visit     = if (has(vis)) sum(is.na(df[[vis]]) | df[[vis]] %in% "") else NA,
      na_date      = if (!is.null(d)) sum(is.na(d)) else NA,
      date_min     = if (!is.null(d) && any(!is.na(d))) as.character(min(d, na.rm = TRUE)) else NA,
      date_max     = if (!is.null(d) && any(!is.na(d))) as.character(max(d, na.rm = TRUE)) else NA,
      phases       = if (!is.null(tt)) paste(names(tt), tt, collapse = "; ") else NA,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, lapply(specs, summ_one))
  show(out); save_csv(out, "01_table_summary.csv")
})

# ---- D. missing codes -----------------------------------------------------
section("D. Missing-value codes (-4 / -5) and NAs in numeric columns", {
  res <- list()
  for (nm in setdiff(names(tb), "DATADIC")) {
    df <- tb[[nm]]
    for (cn in names(df)) {
      x <- df[[cn]]
      if (!is.numeric(x)) next
      n4 <- sum(x == -4, na.rm = TRUE); n5 <- sum(x == -5, na.rm = TRUE); nna <- sum(is.na(x))
      if (n4 + n5 + nna > 0) {
        res[[length(res) + 1]] <- data.frame(table = nm, column = cn, n = length(x),
          n_minus4 = n4, n_minus5 = n5, n_na = nna, stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, res)
  out$pct_missing <- round(100 * (out$n_minus4 + out$n_minus5 + out$n_na) / out$n, 1)
  save_csv(out, "01_missing_codes.csv")
  cat(nrow(out), "numeric columns have NA, -4 or -5. Top 25 by % missing:\n")
  show(out[order(-out$pct_missing), ], 25)
})

# ---- E. diagnosis spine ---------------------------------------------------
dxv <- NULL; ever_mci <- NULL; last_dx <- NULL
section("E. Diagnosis spine (DXSUM)", {
  dx <- tb$DXSUM
  dx$EXAMDATE  <- to_date(dx$EXAMDATE)
  dx$DIAGNOSIS <- suppressWarnings(as.numeric(dx$DIAGNOSIS))
  cat("DIAGNOSIS codes by PHASE (1=CN, 2=MCI, 3=Dementia, 10=TEAM no-diagnosis):\n")
  dx %>% count(PHASE, DIAGNOSIS) %>% show()
  cat("\nRows with missing EXAMDATE:", sum(is.na(dx$EXAMDATE)),
      "| missing DIAGNOSIS:", sum(is.na(dx$DIAGNOSIS)),
      "| DIAGNOSIS == 10:", sum(dx$DIAGNOSIS %in% 10), "\n")
  dxv <- dx %>% filter(!is.na(EXAMDATE), DIAGNOSIS %in% c(1, 2, 3)) %>% arrange(RID, EXAMDATE)
  cat("Valid spine rows:", nrow(dxv), "| subjects:", n_distinct(dxv$RID), "\n")
  cat("Duplicate RID+EXAMDATE:", sum(duplicated(dxv[c("RID", "EXAMDATE")])),
      "| duplicate RID+VISCODE2:", sum(duplicated(dxv[c("RID", "VISCODE2")])), "\n")
  cat("\nVisits per subject (valid diagnosis rows):\n")
  print(summary(as.numeric(table(dxv$RID))))
  cat("\nFirst diagnosis per subject:\n")
  dxv %>% group_by(RID) %>% summarise(first_dx = first(DIAGNOSIS)) %>% count(first_dx) %>% show()
  cat("\nTransitions between consecutive visits (prev -> next):\n")
  tr <- dxv %>% group_by(RID) %>% mutate(prev = lag(DIAGNOSIS)) %>% ungroup() %>%
    filter(!is.na(prev)) %>% count(prev, DIAGNOSIS)
  show(tr); save_csv(tr, "01_dx_transitions.csv")
  cat("\nVISCODE2 frequency (top 25):\n")
  dxv %>% count(VISCODE2, sort = TRUE) %>% show(25)
  ever_mci <- unique(dxv$RID[dxv$DIAGNOSIS == 2])
  last_dx  <- dxv %>% group_by(RID) %>%
    summarise(last_dx_date = max(EXAMDATE), n_visits = n())
  cat("\nSubjects ever diagnosed MCI:", length(ever_mci), "\n")
})

# ---- F. subject-level covariates -----------------------------------------
demo <- NULL
section("F. Subject-level covariates (DM, PTDEMOG, APOERES)", {
  dm <- tb$DM
  dm$RID <- suppressWarnings(as.integer(dm$SUBJID))
  dm <- dm %>% distinct(RID, .keep_all = TRUE)
  pt <- tb$PTDEMOG
  pt$VISDATE <- to_date(pt$VISDATE)
  pt$edu_num <- suppressWarnings(as.numeric(pt$PTEDUCAT))
  edu <- pt %>% arrange(RID, VISDATE) %>% group_by(RID) %>%
    summarise(edu = first(edu_num[!is.na(edu_num)]))
  ap <- tb$APOERES %>% filter(!is.na(GENOTYPE), GENOTYPE != "") %>%
    group_by(RID) %>% summarise(genotype = first(GENOTYPE)) %>%
    mutate(e4 = str_count(genotype, "4"))
  d <- dm %>% transmute(RID, sex = SEX, race_raw = RACE, ethnic = ETHNIC,
                        age_enrol = AGE, enrol_date = to_date(DMDTC),
                        origprot = ORIGPROT, dth = DTHFL %in% "Yes",
                        death_date = to_date(DTHDTC)) %>%
    left_join(edu, by = "RID") %>% left_join(ap, by = "RID")
  d <- d %>% mutate(
    race = case_when(is.na(race_raw) | race_raw == "Unknown" ~ "Unknown",
                     race_raw == "White" ~ "White",
                     race_raw == "Black or African American" ~ "Black",
                     race_raw == "Asian" ~ "Asian",
                     TRUE ~ "Other/multiple"),
    hisp = case_when(ethnic == "Hispanic or Latino" ~ "Hispanic",
                     ethnic == "Not Hispanic or Latino" ~ "Not Hispanic",
                     TRUE ~ "Unknown"),
    edu_grp   = case_when(is.na(edu) ~ "missing", edu <= 12 ~ "<=12y", TRUE ~ ">12y"),
    apoe4_grp = case_when(is.na(e4) ~ "missing", e4 == 0 ~ "0", TRUE ~ "1+"))
  demo <- d
  cat("DM subjects:", nrow(d), "| ever-MCI subjects with a DM row:",
      sum(ever_mci %in% d$RID), "of", length(ever_mci), "\n")
  for (v in c("sex", "race", "hisp", "edu_grp", "apoe4_grp", "origprot")) {
    cat("\n", v, " (ever-MCI subjects):\n", sep = "")
    print(table(d[[v]][d$RID %in% ever_mci], useNA = "ifany"))
  }
  cat("\nDeaths among ever-MCI subjects:", sum(d$dth[d$RID %in% ever_mci]), "\n")
})

# ---- G. alignment of assessments to the diagnosis spine -------------------
section("G. Assessment dates vs diagnosis-visit dates (join by visit code?)", {
  spine <- dxv %>% group_by(RID, VISCODE2) %>% summarise(dx_date = first(EXAMDATE)) %>% ungroup()
  al <- list(list("CDR", "VISDATE"), list("ADAS", "VISDATE"), list("FAQ", "VISDATE"),
             list("MMSE", "VISDATE"), list("NEUROBAT", "VISDATE"),
             list("UCSFFSX7", "EXAMDATE"), list("PLASMA", "EXAMDATE"), list("AV45", "EXAMDATE"))
  out <- lapply(al, function(a) {
    df <- tb[[a[[1]]]]
    if (is.null(df) || !(a[[2]] %in% names(df))) return(NULL)
    x <- data.frame(RID = df$RID, VISCODE2 = df$VISCODE2, d = to_date(df[[a[[2]]]]),
                    stringsAsFactors = FALSE)
    m   <- left_join(x, spine, by = c("RID", "VISCODE2"))
    gap <- as.numeric(m$d - m$dx_date)
    q   <- suppressWarnings(quantile(gap, c(0.05, 0.5, 0.95), na.rm = TRUE))
    data.frame(table = a[[1]], rows = nrow(m),
               matched_by_code = sum(!is.na(m$dx_date)),
               pct_matched = round(100 * mean(!is.na(m$dx_date)), 1),
               within_30d = sum(abs(gap) <= 30, na.rm = TRUE),
               within_90d = sum(abs(gap) <= 90, na.rm = TRUE),
               within_180d = sum(abs(gap) <= 180, na.rm = TRUE),
               beyond_180d = sum(abs(gap) > 180, na.rm = TRUE),
               gap_q05 = unname(q[1]), gap_median = unname(q[2]), gap_q95 = unname(q[3]),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  show(out); save_csv(out, "01_visit_alignment.csv")
  cat("\n(gap = assessment date minus diagnosis date, in days, among rows matched by RID + VISCODE2)\n")
})

# ---- H. enrolment phases and rollovers ------------------------------------
section("H. Enrolment phase structure (ORIGPROT) and rollovers", {
  ph <- dxv %>% group_by(RID) %>% summarise(phases = paste(sort(unique(PHASE)), collapse = "+"))
  ev <- demo %>% select(RID, origprot) %>% inner_join(ph, by = "RID") %>% filter(RID %in% ever_mci)
  tab <- ev %>% count(origprot, phases) %>% arrange(origprot, desc(n))
  show(tab, 100); save_csv(tab, "01_origprot_by_phases_everMCI.csv")
  cat("\nEver-MCI subjects by original enrolment phase:\n")
  print(table(ev$origprot, useNA = "ifany"))
  cat("\nEver-MCI subjects with any ADNI4 or TEAM diagnosis row:",
      sum(grepl("ADNI4|TEAM", ev$phases)), "\n")
})

# ---- I. how follow-up ends ------------------------------------------------
section("I. How follow-up ends: death, discontinuation, withdrawal (ever-MCI)", {
  data_max <- max(dxv$EXAMDATE)
  cat("Latest diagnosis date in the spine:", format(data_max), "\n")
  rg <- tb$REGISTRY
  cat("\nREGISTRY reasons a visit was not done (all rows):\n")
  rg %>% count(RGREASON, sort = TRUE) %>% show(15)
  cat("\nREGISTRY participant status (all rows):\n")
  rg %>% count(PTSTATUS, sort = TRUE) %>% show(10)
  flags <- rg %>% group_by(RID) %>% summarise(
    discontinued = any(grepl("Discontinued", PTSTATUS)),
    withdrew     = any(RGREASON %in% "Withdrawn Consent"))
  ex <- demo %>% filter(RID %in% ever_mci) %>%
    select(RID, origprot, dth, death_date) %>%
    left_join(last_dx, by = "RID") %>% left_join(flags, by = "RID")
  ex$discontinued[is.na(ex$discontinued)] <- FALSE
  ex$withdrew[is.na(ex$withdrew)]         <- FALSE
  ex$exit <- case_when(
    ex$dth                                 ~ "1 death",
    ex$withdrew                            ~ "2 withdrew consent",
    ex$discontinued                        ~ "3 discontinued",
    ex$last_dx_date >= data_max - 365      ~ "4 recent visit (<=12m before data end)",
    TRUE                                   ~ "5 no recorded exit; last visit >12m before data end")
  cat("\nCrude exit classification (priority order shown by the number):\n")
  tt <- ex %>% count(exit) %>% mutate(pct = round(100 * n / sum(n), 1))
  show(tt); save_csv(tt, "01_exit_types.csv")
  tt2 <- ex %>% count(origprot, exit)
  save_csv(tt2, "01_exit_types_by_origprot.csv")
  cat("\nOverlap of flags (death / withdrew / discontinued):\n")
  print(table(death = ex$dth, withdrew = ex$withdrew, discontinued = ex$discontinued))
  gap <- as.numeric(ex$death_date - ex$last_dx_date)[ex$dth]
  cat("\nDays from last diagnosis visit to death (decedents): \n")
  print(summary(gap))
})

# ---- J. landmark cohorts and subgroup counts ------------------------------
all_lm <- NULL
section("J. Landmark cohorts and subgroup counts (CRUDE: no censoring model)", {
  dem1 <- dxv %>% filter(DIAGNOSIS == 3) %>% select(RID, dem_date = EXAMDATE)
  res <- list()
  for (lm in CFG$landmarks) {
    lmk <- dxv %>% filter(VISCODE2 == lm, DIAGNOSIS == 2) %>%
      group_by(RID) %>% slice(1) %>% ungroup() %>% select(RID, lm_date = EXAMDATE)
    nxt <- lmk %>% inner_join(dem1, by = "RID") %>% filter(dem_date > lm_date) %>%
      group_by(RID) %>% summarise(first_dem = min(dem_date))
    d <- lmk %>% left_join(nxt, by = "RID") %>% left_join(last_dx, by = "RID") %>%
      left_join(demo, by = "RID")
    d$landmark <- lm
    res[[lm]]  <- d
  }
  a <- bind_rows(res)
  a$landmark <- factor(a$landmark, levels = CFG$landmarks)
  a$t_event  <- as.numeric(a$first_dem - a$lm_date) / 365.25
  a$fu_years <- as.numeric(a$last_dx_date - a$lm_date) / 365.25
  a$age_lm   <- a$age_enrol + as.numeric(a$lm_date - a$enrol_date) / 365.25
  a$age_grp  <- as.character(cut(a$age_lm, c(-Inf, 65, 75, Inf),
                                 labels = c("<65", "65-74", "75+"), right = FALSE))
  for (h in CFG$horizons_years) {
    a[[paste0("st", h)]] <- ifelse(!is.na(a$t_event) & a$t_event <= h, "event",
                             ifelse(a$fu_years >= h, "followed", "censored"))
  }
  all_lm <- a

  ov <- do.call(rbind, lapply(CFG$landmarks, function(lm) {
    d <- a[a$landmark == lm, ]
    do.call(rbind, lapply(CFG$horizons_years, function(h) {
      st <- d[[paste0("st", h)]]
      data.frame(landmark = lm, horizon_y = h, n_mci = nrow(d),
                 events = sum(st == "event"),
                 followed_no_event = sum(st == "followed"),
                 censored_before_h = sum(st == "censored"),
                 pct_censored = round(100 * mean(st == "censored"), 1))
    }))
  }))
  cat("MCI at each landmark visit (unique subjects), outcome status by horizon:\n")
  cat("  event = dementia diagnosed within h years of the landmark\n")
  cat("  followed_no_event = no dementia, last diagnosis visit at least h years after landmark\n")
  cat("  censored_before_h = no dementia, last visit earlier than h years after landmark\n\n")
  show(ov, 20); save_csv(ov, "01_landmark_overall.csv")

  sg <- function(v, h) {
    d <- a
    d$level <- as.character(d[[v]]); d$level[is.na(d$level)] <- "missing"
    d$st    <- d[[paste0("st", h)]]
    d %>% group_by(landmark, level) %>%
      summarise(n_mci = n(), events = sum(st == "event"),
                followed_no_event = sum(st == "followed"),
                censored_before_h = sum(st == "censored"), .groups = "drop") %>%
      mutate(pct_censored = round(100 * censored_before_h / n_mci, 1),
             subgroup = v, horizon_y = h) %>%
      select(subgroup, level, landmark, horizon_y, everything())
  }
  vars <- c("sex", "age_grp", "edu_grp", "apoe4_grp", "race", "hisp", "origprot")
  subs <- do.call(rbind, lapply(vars, function(v)
    do.call(rbind, lapply(CFG$horizons_years, function(h) sg(v, h)))))
  subs <- as.data.frame(subs) %>% arrange(subgroup, level, horizon_y, landmark)
  save_csv(subs, "01_landmark_subgroups.csv")
  cat("\nSubgroup counts at the 2-year horizon (all horizons are in the CSV):\n")
  show(subs %>% filter(horizon_y == 2), 200)
})

# ---- K. biomarker availability (tier D) -----------------------------------
section("K. Biomarker availability among MCI at each landmark (exact visit-code match)", {
  have <- function(df, lm) unique(df$RID[df$VISCODE2 == lm])
  out <- do.call(rbind, lapply(CFG$landmarks, function(lm) {
    d <- all_lm[all_lm$landmark == lm, ]
    data.frame(landmark = lm, n_mci = nrow(d),
               plasma_at_visit = sum(d$RID %in% have(tb$PLASMA, lm)),
               av45_at_visit   = sum(d$RID %in% have(tb$AV45, lm)),
               freesurfer_at_visit = sum(d$RID %in% have(tb$UCSFFSX7, lm)),
               plasma_any_time = sum(d$RID %in% tb$PLASMA$RID))
  }))
  show(out); save_csv(out, "01_biomarker_availability.csv")
  cat("\nPlasma rows by phase, ever-MCI subjects only:\n")
  tb$PLASMA %>% filter(RID %in% ever_mci) %>% count(PHASE) %>% show()
  cat("\nFreeSurfer OVERALLQC (all rows):\n")
  print(table(tb$UCSFFSX7$OVERALLQC, useNA = "ifany"))
})

# ---- L. variable check and data dictionary --------------------------------
section("L. Variables needed for the visit table, and DATADIC", {
  dd <- tb$DATADIC
  cat("DATADIC rows:", nrow(dd), "| tables:", n_distinct(dd$TBLNAME), "\n")
  chk <- function(tbl, cols) {
    df <- tb[[tbl]]
    data.frame(table = tbl, column = cols,
               present = if (is.null(df)) NA else cols %in% names(df))
  }
  want <- rbind(
    chk("UCSFFSX7", c("ST29SV", "ST88SV", "ST24CV", "ST83CV", "ST26CV", "ST85CV",
                      "ST37SV", "ST96SV", "ST10CV", "OVERALLQC", "HIPPOQC", "VENTQC", "FSVER")),
    chk("NEUROBAT", c("AVTOT1", "AVTOT2", "AVTOT3", "AVTOT4", "AVTOT5", "AVDEL30MIN", "AVDELTOT")),
    chk("CDR", c("CDRSB", "CDGLOBAL")),
    chk("ADAS", c("TOTSCORE", "TOTAL13")),
    chk("FAQ", "FAQTOTAL"),
    chk("MMSE", "MMSCORE"),
    chk("PTDEMOG", c("PTEDUCAT", "PTGENDER", "PTRACCAT", "PTETHCAT")),
    chk("DM", c("AGE", "SEX", "RACE", "ETHNIC", "DTHFL", "DTHDTC")))
  show(want, 100); save_csv(want, "01_variable_check.csv")
  cat("\nMMSE columns:\n"); print(names(tb$MMSE))
  cat("\nDXSUM DIAGNOSIS code definitions in DATADIC:\n")
  dd %>% filter(TBLNAME == "DXSUM", FLDNAME == "DIAGNOSIS") %>% distinct(PHASE, CODE) %>% show()
  spec <- dd %>% filter(TBLNAME == "ADNIMERGE") %>% distinct(FLDNAME, TEXT, CODE, UNITS)
  cat("\nADNIMERGE rows in DATADIC (the merged-table spec):", nrow(spec), "\n")
  save_csv(spec, "01_adnimerge_spec.csv")
  show(spec %>% filter(str_detect(FLDNAME,
    "RAVLT|ADAS|CDRSB|FAQ|MMSE|PACC|ICV|Hippo|Entorh|Ventric|WholeBrain|Fusiform|MidTemp|APOE4|^DX|^AGE|^PT")), 80)
})

cat("\n", strrep("=", 72), "\nSession info\n", strrep("=", 72), "\n", sep = "")
print(sessionInfo())
sink()
cat("\nDone. Report: ", report, "\nTables:  ", file.path(CFG$dir_out, "audit"), "\n", sep = "")
