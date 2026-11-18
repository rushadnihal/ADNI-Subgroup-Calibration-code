# ADNI subgroup calibration: analysis code

R code for a landmark audit of MCI-to-dementia risk prediction in the Alzheimer's Disease Neuroimaging Initiative (ADNI): calibration by landmark time, protocol era and demographic subgroup, how many events each subgroup contributes, differential loss to follow-up, and record-level versus subject-level splitting.

**Status:** working draft; the accompanying manuscript has not been submitted or peer reviewed.

## Preliminary findings (unpublished)

> **These results are unpublished and have not been peer reviewed.** They come from a working draft, from a single ADNI download (11 November 2026), and may change. Do not cite them as established findings. Only aggregate summaries are shown; no participant-level data are included.

The audit follows 1,194 people with MCI across four landmarks (baseline, 12, 24 and 36 months after baseline), giving 2,732 subject-landmark rows. At baseline there are 1,147 rows and 201 dementia events within 2 years. Out-of-fold predictions use subject-level cross-validation. O/E is observed over expected events, so values below 1 mean the model over-predicts risk.

1. **Overall, the models look well calibrated.** A Cox model with 8 terms (Tier B) has a baseline 2-year AUC of 0.85 (95% CI 0.83 to 0.88) and O/E of 1.01 (0.90 to 1.13). A confirmed-dementia-only outcome gives the same picture (AUC 0.86, O/E 0.97).
2. **Calibration differs by protocol era even when discrimination does not.** For the same Tier B model at baseline, 2-year O/E is 1.34 (1.17 to 1.53) in ADNI1, 0.85 (0.67 to 1.05) in ADNI-GO/2 and 0.50 (0.27 to 0.73) in ADNI3, so it under-predicts early and over-predicts late. A model trained on ADNI1 only and applied to ADNI-GO, ADNI2 and ADNI3 had an unchanged AUC (about 0.88) but a worse O/E (0.53, 0.43 to 0.64, versus 0.74 for the pooled model).
3. **Loss to follow-up is large and could explain some of the gaps.** Of people observable at the 2-year horizon, 17% were lost at the baseline landmark and 35% at the 36-month landmark. Bounds on the missing outcomes show that censoring could reconcile the gaps in most groups (ADNI1 is already above predicted). For the 3-term baseline model (age, APOE4, CDR-SB), a covariate-based censoring correction (IPCW) leaves the era O/E values nearly unchanged, so measured covariates do not explain them.
4. **Most subgroups are too small to audit.** Detecting a 20% miscalibration needs roughly 135 to 184 expected events in a subgroup. None of the 74 subgroup-by-landmark levels examined reaches that. Black and Hispanic participants have fewer than 10 events each at baseline, so only very large miscalibration (O/E below about 0.13 or above about 2.3) could be detected.
5. **Richer models add little.** A 12-term Cox model, elastic-net Cox and a random survival forest improve baseline 2-year AUC by about 0.01 over the Tier B Cox model (the interval excludes 0 for the Cox and elastic-net models, not for the forest). At the later landmarks every interval includes 0.
6. **Record-level splitting is a risk for flexible models only.** Stacking landmark rows and splitting by record leaves about 80% of test rows with a same-person row in training. This inflates AUC by at most 0.002 for the Cox model (interval includes 0) but by 0.014 to 0.022 for random survival forests.

Details for each point are in the scripts (`R/`, listed in `R/README.md`) and the results they write.

## What is and is not in this repository

- **Included:** the analysis code only (`R/`).
- **Not included:** ADNI data, any participant-level or derived participant-level file, and the full result tables and figures (only the summary below is shown). Under the ADNI Data Use Agreement, redistribution of participant-level data is prohibited.

## To reproduce

1. Apply for ADNI access at adni.loni.usc.edu and download the ADNIMERGE2 tables.
2. Place them in a folder called `data and methods/` inside the project folder, with the scripts in `R/`.
3. Install the pinned package versions from [`renv.lock`](renv.lock) (R 4.6.2): `Rscript -e 'install.packages("renv"); renv::activate(); renv::restore(prompt = FALSE)'`. CI installs the same versions and checks that every script parses; it cannot run the analysis without ADNI data.
4. From the project folder run `source("R/run_all.R")` (about 25 minutes). See `R/README.md` for what each script does and writes.

The scripts were developed on data downloaded on 11 November 2026 (diagnosis visits to 2026-11-09). ADNI data are updated over time, so later downloads may give slightly different numbers. Each script has built-in checks that print PASS or FAIL, and `run_all.R` summarises them.

## ADNI acknowledgement

Data used in the preparation of this article were obtained from the Alzheimer's Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). The ADNI was launched in 2003 as a public-private partnership, led by Principal Investigator Michael W. Weiner, MD. The original goal of ADNI was to test whether serial magnetic resonance imaging (MRI), positron emission tomography (PET), other biological markers, and clinical and neuropsychological assessment can be combined to measure the progression of mild cognitive impairment (MCI) and early Alzheimer's disease (AD). The current goals include validating biomarkers for clinical trials, improving the generalizability of ADNI data by increasing diversity in the participant cohort, and to provide data concerning the diagnosis and progression of Alzheimer’s disease to the scientific community. For up-to-date information, see adni.loni.usc.edu.

Data collection and sharing for the Alzheimer's Disease Neuroimaging Initiative (ADNI) is funded by the National Institute on Aging (National Institutes of Health Grant U19AG024904). The grantee organization is the Northern California Institute for Research and Education. In the past, ADNI has also received funding from the National Institute of Biomedical Imaging and Bioengineering, the Canadian Institutes of Health Research, and private sector contributions through the Foundation for the National Institutes of Health (FNIH) including generous contributions from the following: AbbVie, Alzheimer’s Association; Alzheimer’s Drug Discovery Foundation; Araclon Biotech; BioClinica, Inc.; Biogen; Bristol- Myers Squibb Company; CereSpir, Inc.; Cogstate; Eisai Inc.; Elan Pharmaceuticals, Inc.; Eli Lilly and Company; EuroImmun; F. Hoffmann-La Roche Ltd and its affiliated company Genentech, Inc.; Fujirebio; GE Healthcare; IXICO Ltd.; Janssen Alzheimer Immunotherapy Research & Development, LLC.; Johnson & Johnson Pharmaceutical Research & Development LLC.; Lumosity; Lundbeck; Merck & Co., Inc.; Meso Scale Diagnostics, LLC.; NeuroRx Research; Neurotrack Technologies; Novartis Pharmaceuticals Corporation; Pfizer Inc.; Piramal Imaging; Servier; Takeda Pharmaceutical Company; and Transition Therapeutics.

## Licence

Code: MIT (see `LICENSE`).
