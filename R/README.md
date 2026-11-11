# R scripts

Run order. Each script sources `00_setup.R` itself, so you only ever run the numbered one.

| Script | Purpose | Writes to |
| --- | --- | --- |
| `00_setup.R` | Paths, seed, package check, helper functions. Not run on its own. | `output/` folders |
| `01_audit.R` | Read-only data audit. No modelling. | `output/audit/` |
| `02_build_visits.R` | One row per diagnosis visit, with the latest assessment on or before each visit date. Leakage guard must print PASS. | `output/interim/visits.rds`, `output/audit/02_*` |
| `03_outcomes.R` | Landmark cohort (MCI at bl, m12, m24, m36), dementia and death outcomes, censoring by subgroup. | `output/interim/outcomes.rds`, `output/audit/03_*` |
| `04_splits_baseline.R` | Subject-level 5-fold split, leakage and synthetic tests (each prints PASS or FAIL), cross-validated baseline Cox, calibration by subgroup. | `output/interim/folds.rds`, `output/audit/04_*` |
| `05_tiers_bootstrap.R` | Feature tiers A, B, C (plus protocol-era sensitivity), cross-validated on the same folds, with subject-level bootstrap intervals for AUC, Brier, O/E and subgroup O/E. About 4 to 8 minutes. | `output/audit/05_*` |
| `06_censoring_power.R` | Could loss to follow-up explain the ADNI3 and race calibration results? Risk bounds with tipping fraction, covariate-adjusted IPCW, and a power table for subgroup calibration. About 1 to 2 minutes. | `output/audit/06_*` |
| `07_model_classes.R` | Benchmark of model classes on the same folds, features and metrics: Cox (tiers B and C), elastic-net Cox, random survival forest, each with and without a protocol-era term. Bootstrap intervals, AUC gain over Cox-B, subgroup O/E. Needs `glmnet` and `ranger`. About 8 to 15 minutes. | `output/interim/oof_predictions.rds`, `output/audit/07_*` |
| `08_calibration_over_time.R` | Subgroup O/E at every landmark, horizon and model from the saved predictions, bootstrap intervals, ratios of O/E between subgroup levels, and two figures. About 2 to 4 minutes. | `output/audit/08_*`, `output/figures/08_*.png` |
| `09_decision_curves_conformal.R` | Decision curves (net benefit by risk threshold, whole cohort and by subgroup) and cross-conformal prediction sets with IPCW coverage by subgroup, from the saved predictions. A few minutes. | `output/audit/09_*`, `output/figures/09_*.png` |
| `10_plasma_competing_risk.R` | Tier D: who has plasma, Tier B with and without plasma on the same rows (bl, m12), shuffled-plasma control; death as a competing event (Aalen-Johansen vs 1 - KM, cause-specific Cox, Fine-Gray). A few minutes. | `output/audit/10_*`, `output/figures/10_*.png` |
| `11_events_needed.R` | Events needed to detect a 20% O/E deviation (and group contrasts) vs events ADNI has per subgroup and landmark. Under a minute. | `output/audit/11_*`, `output/figures/11_*.png` |
| `12_paper_figures.R` | Paper-ready figures (PNG and PDF) redrawn from the 08 to 11 CSVs: clipped intervals marked, decision-curve gain with intervals, paired plasma AUC differences, subgroup O/E over the undetectable range. Under a minute. | `output/figures/12_*`, `output/audit/12_*` |
| `13_sensitivity.R` | Sensitivity analyses: (A) confirmed-dementia outcome (metrics, subgroup and era O/E, three scenarios, paired bootstrap); (B) censoring bounds, tipping fractions and IPCW check repeated for Cox-B. A few minutes. | `output/audit/13_*`, `output/figures/13_*` |
| `14_temporal_holdout.R` | (A) Cox-B trained on ADNI1 only and applied to ADNI-GO, ADNI2 and ADNI3 (O/E and AUC against the pooled out-of-fold predictions, paired bootstrap); (B) record-level versus subject-level splits on stacked landmark rows (AUC inflation). About 15 to 25 minutes. | `output/audit/14_*`, `output/figures/14_*` |
| `15_forest_split_check.R` | Record-level versus subject-level splits on stacked landmark rows for a random survival forest (default and deep), compared with Cox-B from 14. About 5 to 15 minutes. | `output/audit/15_*`, `output/figures/15_*` |
| `run_all.R` | Runs 01 to 15 in order, each in a fresh R process, and writes a summary of exit status, minutes, checks passed and failed. About 25 minutes on the author's machine (23 minutes measured; allow up to an hour). `source("R/run_all.R")` or `Rscript R/run_all.R`. | `output/audit/run_all_summary.csv`, `output/logs/run_NN.log`, `output/logs/session_info.txt` |

`utils_eval.R` holds shared helper functions for 06 onward. It is sourced automatically and is not run on its own.

## How to run

1. Set the working directory to the project folder (`analysis1111`), or to its `R/` subfolder.
2. Start a fresh R session, then `source("R/01_audit.R")`, then 02 to 15 (each reads the previous one's output).
3. If it stops with "Install missing packages", run the `install.packages(...)` line it prints, then repeat.

Needs R 4.x with `dplyr` (1.0 or later), `readr` and `stringr`. Run time is about 1 to 3 minutes.

## What comes back

`01_audit.R` writes a full text log, `output/audit/01_audit_report.txt`, plus small CSV tables. They hold counts only, with no subject IDs. If a section fails, the log says `!! SECTION FAILED` with the reason and the rest still runs.

## Reference counts to check the R output against

These were computed independently in Python from the same files on 11 Nov 2026.

- Diagnosis table: 16,640 rows, 110 with no exam date, 46 with no diagnosis, 94 with diagnosis code 10.
- Valid spine: 16,437 rows, 3,685 subjects, 37 duplicate subject-and-date pairs, 2 duplicate subject-and-visit-code pairs.
- Subjects ever diagnosed MCI: 1,759.
- MCI at baseline: 1,287 subjects. Dementia within 1, 2 and 3 years: 80, 202 and 278. Censored before the horizon: 215, 383 and 476.
- MCI at month 36: 449 subjects. Dementia within 1, 2 and 3 years: 16, 41 and 51.

Small differences (a handful of subjects) can come from how duplicates are handled.

## Rules

- Do not copy anything from `data and methods/` or `output/` into a repository (see `.gitignore`).
- Every script writes only inside `output/`. Nothing in `data and methods/` is modified.
