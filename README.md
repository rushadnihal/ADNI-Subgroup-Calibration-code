# ADNI subgroup calibration: analysis code

R code for a landmark audit of MCI-to-dementia risk prediction in the Alzheimer's Disease Neuroimaging Initiative (ADNI): calibration by landmark time, protocol era and demographic subgroup, how many events each subgroup contributes, differential loss to follow-up, and record-level versus subject-level splitting.

**Status:** working draft; the accompanying manuscript has not been submitted or peer reviewed.

## What is and is not in this repository

- **Included:** the analysis code only (`R/`).
- **Not included:** ADNI data, any participant-level or derived participant-level file, and results. Under the ADNI Data Use Agreement, redistribution of participant-level data is prohibited.

## To reproduce

1. Apply for ADNI access at adni.loni.usc.edu and download the ADNIMERGE2 tables.
2. Place them in a folder called `data and methods/` inside the project folder, with the scripts in `R/`.
3. From the project folder run `source("R/run_all.R")` (about 25 minutes). See `R/README.md` for what each script does and writes.

The scripts were developed on data downloaded on 11 November 2026 (diagnosis visits to 2026-11-09). ADNI data are updated over time, so later downloads may give slightly different numbers. Each script has built-in checks that print PASS or FAIL, and `run_all.R` summarises them.

## ADNI acknowledgement

Data used in the preparation of this article were obtained from the Alzheimer's Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). The ADNI was launched in 2003 as a public-private partnership, led by Principal Investigator Michael W. Weiner, MD. The original goal of ADNI was to test whether serial magnetic resonance imaging (MRI), positron emission tomography (PET), other biological markers, and clinical and neuropsychological assessment can be combined to measure the progression of mild cognitive impairment (MCI) and early Alzheimer's disease (AD). The current goals include validating biomarkers for clinical trials, improving the generalizability of ADNI data by increasing diversity in the participant cohort, and to provide data concerning the diagnosis and progression of Alzheimer’s disease to the scientific community. For up-to-date information, see adni.loni.usc.edu.

Data collection and sharing for the Alzheimer's Disease Neuroimaging Initiative (ADNI) is funded by the National Institute on Aging (National Institutes of Health Grant U19AG024904). The grantee organization is the Northern California Institute for Research and Education. In the past, ADNI has also received funding from the National Institute of Biomedical Imaging and Bioengineering, the Canadian Institutes of Health Research, and private sector contributions through the Foundation for the National Institutes of Health (FNIH) including generous contributions from the following: AbbVie, Alzheimer’s Association; Alzheimer’s Drug Discovery Foundation; Araclon Biotech; BioClinica, Inc.; Biogen; Bristol- Myers Squibb Company; CereSpir, Inc.; Cogstate; Eisai Inc.; Elan Pharmaceuticals, Inc.; Eli Lilly and Company; EuroImmun; F. Hoffmann-La Roche Ltd and its affiliated company Genentech, Inc.; Fujirebio; GE Healthcare; IXICO Ltd.; Janssen Alzheimer Immunotherapy Research & Development, LLC.; Johnson & Johnson Pharmaceutical Research & Development LLC.; Lumosity; Lundbeck; Merck & Co., Inc.; Meso Scale Diagnostics, LLC.; NeuroRx Research; Neurotrack Technologies; Novartis Pharmaceuticals Corporation; Pfizer Inc.; Piramal Imaging; Servier; Takeda Pharmaceutical Company; and Transition Therapeutics.

## Licence

Code: MIT (see `LICENSE`).
