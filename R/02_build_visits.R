# 02_build_visits.R --------------------------------------------------------
# Builds the visit-level table: one row per valid diagnosis visit, with
# covariates and the MOST RECENT assessment recorded ON OR BEFORE that visit's
# date ("as-of" join). Nothing recorded after a visit can enter that visit's
# row, which is what the outline's "inputs only at or before the landmark"
# rule requires.
#
# Writes: output/interim/visits.rds            (subject-level, stays local)
#         output/audit/02_*.csv, 02_build_visits_report.txt  (counts only)
# Run:    source("R/02_build_visits.R")   (needs dplyr >= 1.1)
# -------------------------------------------------------------------------

if (!exists("CFG")) source(if (file.exists("R/00_setup.R")) "R/00_setup.R" else "00_setup.R")

report <- file.path(CFG$dir_out, "audit", "02_build_visits_report.txt")
sink(report, split = TRUE)
options(width = 170)
cat("Visit table build | run:", format(Sys.time()), "\n")

# ---- settings you may want to change --------------------------------------
# Maximum age (days) of an assessment that may be carried into a visit.
# Defaults come from the audit: scales are usually recorded within ~100 days
# before the diagnosis date; imaging and plasma are sparser.
# Widened after the first run: baseline CDR-SB is often recorded well before
# the diagnosis date (82% within 90 d, 96% within 180 d), and MRI is sparse at
# later landmarks (m36: 35% within 90 d, 90% within 730 d).
LOOKBACK <- c(cdr = 180, adas = 90, faq = 90, mmse = 180, ravlt = 90,
              fs = 730, plasma = 365)
# ---------------------------------------------------------------------------

# numeric conversion that turns ADNI missing codes (-4, -5) into NA
num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% CFG$missing_codes] <- NA
  x
}
step <- function(label, n) cat(sprintf("  %-66s %7d\n", label, n))

# ---- 1. diagnosis spine ---------------------------------------------------
spine <- NULL
section("1. Diagnosis spine and exclusion log", {
  dx <- read_adni_csv("DXSUM_11Nov2026.csv")
  dx$EXAMDATE  <- to_date(dx$EXAMDATE)
  dx$DIAGNOSIS <- suppressWarnings(as.numeric(dx$DIAGNOSIS))
  step("input rows", nrow(dx))
  step("excluded: missing EXAMDATE", sum(is.na(dx$EXAMDATE)))
  step("excluded: DIAGNOSIS missing or code 10 (TEAM, no diagnosis)",
       sum(!is.na(dx$EXAMDATE) & !(dx$DIAGNOSIS %in% c(1, 2, 3))))
  v <- dx %>% filter(!is.na(EXAMDATE), DIAGNOSIS %in% c(1, 2, 3))
  conflicts <- v %>% group_by(RID, EXAMDATE) %>%
    summarise(n = n(), nd = n_distinct(DIAGNOSIS)) %>% ungroup() %>% filter(nd > 1)
  step("excluded: same subject+date with CONFLICTING diagnoses (rows)",
       nrow(semi_join(v, conflicts, by = c("RID", "EXAMDATE"))))
  v <- anti_join(v, conflicts, by = c("RID", "EXAMDATE"))
  n_before <- nrow(v)
  v <- v %>% arrange(RID, EXAMDATE) %>% group_by(RID, EXAMDATE) %>% slice(1) %>% ungroup()
  step("excluded: exact duplicate subject+date rows", n_before - nrow(v))
  spine <- v %>% select(RID, PTID, PHASE, VISCODE2, EXAMDATE, DIAGNOSIS)
  step("ROWS IN SPINE", nrow(spine))
  step("SUBJECTS IN SPINE", n_distinct(spine$RID))
})

# ---- 2. subject-level covariates -----------------------------------------
covs <- NULL
section("2. Subject-level covariates", {
  dm <- load_rda("DM")
  dm$RID <- suppressWarnings(as.integer(dm$SUBJID))
  dm <- dm %>% distinct(RID, .keep_all = TRUE)
  pt <- load_rda("PTDEMOG")
  pt$VISDATE <- to_date(pt$VISDATE)
  pt$edu_num <- num(pt$PTEDUCAT)
  edu <- pt %>% arrange(RID, VISDATE) %>% group_by(RID) %>%
    summarise(edu = first(edu_num[!is.na(edu_num)]))
  ap <- load_rda("APOERES") %>% filter(!is.na(GENOTYPE), GENOTYPE != "") %>%
    group_by(RID) %>% summarise(apoe_genotype = first(GENOTYPE)) %>%
    mutate(apoe4 = str_count(apoe_genotype, "4"))
  covs <- dm %>% transmute(RID, sex = SEX, race_raw = RACE, ethnic = ETHNIC,
                           age_enrol = AGE, enrol_date = to_date(DMDTC),
                           origprot = ORIGPROT, died = DTHFL %in% "Yes",
                           death_date = to_date(DTHDTC)) %>%
    left_join(edu, by = "RID") %>% left_join(ap, by = "RID") %>%
    mutate(
      race = case_when(is.na(race_raw) | race_raw == "Unknown" ~ "Unknown",
                       race_raw == "White" ~ "White",
                       race_raw == "Black or African American" ~ "Black",
                       race_raw == "Asian" ~ "Asian",
                       TRUE ~ "Other/multiple"),
      hisp = case_when(ethnic == "Hispanic or Latino" ~ "Hispanic",
                       ethnic == "Not Hispanic or Latino" ~ "Not Hispanic",
                       TRUE ~ "Unknown"),
      edu_grp   = case_when(is.na(edu) ~ "missing", edu <= 12 ~ "<=12y", TRUE ~ ">12y"),
      apoe4_grp = case_when(is.na(apoe4) ~ "missing", apoe4 == 0 ~ "0", TRUE ~ "1+"))
  cat("Subjects with covariates:", nrow(covs), "\n")
  cat("Spine subjects missing from DM:", sum(!(spine$RID %in% covs$RID)), "\n")
})

# ---- 3. assessment sources (long format: RID, date, features) -------------
src <- list()
finish_src <- function(df) {
  feats <- setdiff(names(df), c("RID", "date"))
  keep  <- !is.na(df$RID) & !is.na(df$date) & rowSums(!is.na(df[feats])) > 0
  df[keep, ] %>% arrange(RID, date)
}
section("3. Assessment sources", {
  x <- load_rda("CDR")
  src$cdr <- finish_src(data.frame(RID = x$RID, date = to_date(x$VISDATE),
                                   cdrsb = num(x$CDRSB), cdr_global = num(x$CDGLOBAL)))
  x <- load_rda("ADAS")
  src$adas <- finish_src(data.frame(RID = x$RID, date = to_date(x$VISDATE),
                                    adas13 = num(x$TOTAL13), adas11 = num(x$TOTSCORE)))
  x <- load_rda("FAQ")
  src$faq <- finish_src(data.frame(RID = x$RID, date = to_date(x$VISDATE),
                                   faq = num(x$FAQTOTAL)))
  x <- load_rda("MMSE")
  src$mmse <- finish_src(data.frame(RID = x$RID, date = to_date(x$VISDATE),
                                    mmse = num(x$MMSCORE)))
  x <- read_adni_csv("NEUROBAT_11Nov2026.csv")
  src$ravlt <- finish_src(data.frame(
    RID = x$RID, date = to_date(x$VISDATE),
    ravlt_immediate  = num(x$AVTOT1) + num(x$AVTOT2) + num(x$AVTOT3) + num(x$AVTOT4) + num(x$AVTOT5),
    ravlt_learning   = num(x$AVTOT5) - num(x$AVTOT1),
    ravlt_forgetting = num(x$AVTOT5) - num(x$AVDEL30MIN)))

  fs <- read_adni_csv("UCSFFSX7_11Nov2026.csv")
  fs$date <- to_date(fs$EXAMDATE)
  fs$run  <- to_date(fs$RUNDATE)
  step("FreeSurfer rows in", nrow(fs))
  step("  excluded: OVERALLQC == Fail", sum(fs$OVERALLQC %in% "Fail"))
  fs <- fs[!(fs$OVERALLQC %in% "Fail"), ]
  n_before <- nrow(fs)
  fs <- fs %>% arrange(RID, date, desc(run)) %>% group_by(RID, date) %>% slice(1) %>% ungroup()
  step("  excluded: repeat processing of the same scan date (kept latest run)", n_before - nrow(fs))
  src$fs <- finish_src(data.frame(
    RID = fs$RID, date = fs$date,
    hippocampus = num(fs$ST29SV) + num(fs$ST88SV),
    entorhinal  = num(fs$ST24CV) + num(fs$ST83CV),
    fusiform    = num(fs$ST26CV) + num(fs$ST85CV),
    mid_temporal = num(fs$ST40CV) + num(fs$ST99CV),
    ventricles  = num(fs$ST37SV) + num(fs$ST96SV),
    icv         = num(fs$ST10CV),
    fs_version  = as.character(fs$FSVER), stringsAsFactors = FALSE))

  pl <- read_adni_csv("UPENN_PLASMA_FUJIREBIO_QUANTERIX_11Nov2026.csv")
  pl$date <- to_date(pl$EXAMDATE)
  n_before <- nrow(pl)
  pl <- pl %>% arrange(RID, date) %>% group_by(RID, date) %>% slice(1) %>% ungroup()
  step("Plasma rows in", n_before)
  step("  excluded: repeat rows for the same subject+date", n_before - nrow(pl))
  src$plasma <- finish_src(data.frame(
    RID = pl$RID, date = pl$date,
    ptau217 = num(pl$pT217_F), ab42_ab40 = num(pl$AB42_AB40_F),
    nfl = num(pl$NfL_Q), gfap = num(pl$GFAP_Q)))

  cat("\nSource rows after cleaning:\n")
  for (nm in names(src)) step(nm, nrow(src[[nm]]))
})

# ---- 4. how far back must we look? (informs LOOKBACK) ----------------------
section("4. Coverage at MCI landmark visits by look-back window", {
  lm <- spine %>% filter(VISCODE2 %in% CFG$landmarks, DIAGNOSIS == 2) %>%
    group_by(RID, VISCODE2) %>% slice(1) %>% ungroup() %>%
    transmute(RID, landmark = VISCODE2, lm_date = EXAMDATE)
  gaps <- function(s) {
    j <- lm %>% inner_join(s[c("RID", "date")], by = "RID", relationship = "many-to-many") %>%
      filter(date <= lm_date) %>% mutate(gap = as.numeric(lm_date - date)) %>%
      group_by(RID, landmark) %>% summarise(gap = min(gap)) %>% ungroup()
    lm %>% select(RID, landmark) %>% left_join(j, by = c("RID", "landmark")) %>% pull(gap)
  }
  ws  <- c(0, 30, 90, 180, 365, 730)
  out <- do.call(rbind, lapply(names(src), function(nm) {
    g <- gaps(src[[nm]])
    do.call(rbind, lapply(CFG$landmarks, function(l) {
      gg <- g[lm$landmark == l]
      r  <- data.frame(source = nm, landmark = l, n_mci = length(gg))
      for (w in ws) r[[paste0("within_", w, "d")]] <- sum(gg <= w, na.rm = TRUE)
      r$any_prior <- sum(!is.na(gg))
      r
    }))
  }))
  show(out, 100); save_csv(out, "02_asof_coverage.csv")
  cat("\n(counts of MCI landmark visits that have an assessment recorded on or before the visit date, within the window)\n")
})

# ---- 5. as-of join onto every spine row ------------------------------------
asof <- function(sp, s, lookback) {
  feats <- setdiff(names(s), c("RID", "date"))
  j <- sp[c("RID", "EXAMDATE")] %>%
    inner_join(s, by = "RID", relationship = "many-to-many") %>%
    filter(date <= EXAMDATE, as.numeric(EXAMDATE - date) <= lookback) %>%
    arrange(RID, EXAMDATE, desc(date)) %>%
    group_by(RID, EXAMDATE) %>% slice(1) %>% ungroup()
  j$age_days <- as.numeric(j$EXAMDATE - j$date)
  names(j)[names(j) == "date"]     <- "src_date"
  names(j)[names(j) == "age_days"] <- "age_days"
  list(tab = j, feats = feats)
}

visits <- NULL
section("5. Assemble the visit table", {
  visits <- spine %>% left_join(covs, by = "RID")
  visits$age_visit <- visits$age_enrol + as.numeric(visits$EXAMDATE - visits$enrol_date) / 365.25
  guard <- list()
  for (nm in names(src)) {
    a <- asof(spine, src[[nm]], LOOKBACK[[nm]])
    t <- a$tab
    keep <- c("RID", "EXAMDATE", a$feats, "src_date", "age_days")
    t <- t[keep]
    names(t)[names(t) == "src_date"] <- paste0(nm, "_date")
    names(t)[names(t) == "age_days"] <- paste0(nm, "_age_days")
    visits <- left_join(visits, t, by = c("RID", "EXAMDATE"))
    guard[[nm]] <- min(t[[paste0(nm, "_age_days")]], na.rm = TRUE)
    cat(sprintf("  %-8s look-back %3d d | matched %6d of %6d visits (%.1f%%)\n",
                nm, LOOKBACK[[nm]], nrow(t), nrow(spine), 100 * nrow(t) / nrow(spine)))
  }
  cat("\nLeakage guard: smallest (visit date - assessment date) per source, must be >= 0:\n")
  print(unlist(guard))
  cat(if (all(unlist(guard) >= 0)) "PASS: no assessment dated after its visit.\n" else "FAIL: an assessment is dated after its visit!\n")
  cat("\nRows:", nrow(visits), "| columns:", ncol(visits), "| subjects:", n_distinct(visits$RID), "\n")
  saveRDS(visits, file.path(CFG$dir_out, "interim", "visits.rds"))
  cat("Saved output/interim/visits.rds\n")
})

# ---- 6. feature coverage ---------------------------------------------------
section("6. Feature coverage: all visits vs MCI landmark visits", {
  feats <- unlist(lapply(src, function(s) setdiff(names(s), c("RID", "date", "fs_version"))))
  feats <- c("age_visit", "edu", "apoe4", feats)
  is_lm <- visits$VISCODE2 %in% CFG$landmarks & visits$DIAGNOSIS == 2
  out <- data.frame(feature = feats,
    pct_nonmissing_all_visits = sapply(feats, function(f) round(100 * mean(!is.na(visits[[f]])), 1)),
    n_mci_landmark_visits = sum(is_lm),
    pct_nonmissing_mci_landmarks = sapply(feats, function(f) round(100 * mean(!is.na(visits[[f]][is_lm])), 1)),
    row.names = NULL)
  show(out, 60); save_csv(out, "02_feature_coverage.csv")
})

# ---- 7. fields still to map ------------------------------------------------
section("7. Candidate FreeSurfer fields for middle temporal and whole brain (from DATADIC)", {
  dd <- read_adni_csv("DATADIC_11Nov2026.csv")
  cand <- dd %>% filter(TBLNAME == "UCSFFSX7",
                        str_detect(TEXT, regex("MiddleTemporal|Middle Temporal|BrainSeg|WholeBrain|Whole Brain|Brain Segmentation", ignore_case = TRUE))) %>%
    distinct(FLDNAME, TEXT)
  show(cand, 40); save_csv(cand, "02_fs_field_candidates.csv")
})

sink()
cat("\nDone. Report: ", report, "\n", sep = "")
