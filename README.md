# Large-Scale, Multi-Ancestry GWAS Modeling for Parkinson's Disease Motor Progression using Hoehn & Yahr Staging

**Last updated:** October 19, 2026
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24204282.svg)](https://doi.org/10.5281/zenodo.24204282)
---

## Overview
This repository accompanies the manuscript entitled: "**Large-Scale, Multi-Ancestry GWAS Modeling for Parkinson's Disease Motor Progression using Hoehn & Yahr Staging**". 

The study serves as a framework for modeling multi-ancestry, multi-cohort clinical data for motor progression. This framework using Hoehn & Yahr data establishes modeling for future progression GWAS. 

---

### Summary
- HY staging is a commonly used clinical scale to assess the severity of Parkinson's disease, ranging from stage 1 (mild) to stage 5 (severe). Reaching HY stage 3 (HY3) is a widely accepted milestone for PD.
- This pipeline tests HY3 as a stable outcome, then estimates how age, sex and disease duration relate to HY3 across GP2 cohorts. It uses two models:
    - **Logistic regression (LR)** on cross-sectional data, reported as odds ratios.
    - **Cox proportional hazards (CPH)** on longitudinal survival data, reported as hazard ratios.
- Each model is fit within each cohort, and the cohort estimates are combined with a random-effects meta-analysis. The pipeline has the following scripts and paths are read from `config.yaml`.:
1. `01_ppmi_create_datasets.py` creates dataset for HY3 validation.
2. `02_hy_ppmi_analysis.R` runs survival analysis for ON/OFF clinics.

    - **ON vs OFF visits:** HY stage during ON versus OFF assessment in the same visit
    - **On Clinic:** HY assessed during ON state; falls back to Unknown if ON unavailable
    - **Off Clinic:** HY assessed during OFF state; falls back to Unknown if OFF unavailable
    - **Logistic regression:** HY ≥ 3 at baseline ~ age + sex + disease duration
    - **Survival analysis:** Time to HY ≥ 3, Cox PH model

3. `03_make_dataset.py` prepares and describes the GP2 HY data.

    Combines GP2 cross-sectional and longitudinal HY records into one harmonized dataset:

    - Removes duplicate visits and drops cohorts with fewer than 20 participants.
    - Derives disease duration at each visit, age at onset (raw and mean-centered) and sex labels.
    - Labels each cohort as longitudinal or cross-sectional (longitudinal if any participant has more than one visit).
    - Produces two versions of the data:
    - The full dataset.
    - A cleaned analysis dataset, which:
        - excludes disease duration over 20 years and HY = 0
        - excludes HY ≥ 4 within the first year, which suggests atypical parkinsonism
        - collapses half-stages (1.5→1, 2.5→2, and so on)
    - Writes descriptive figures and a per-cohort baseline summary table for each version:
    - disease duration by HY stage
    - HY versus duration by sex
    - age at onset by cohort
    - sample sizes
    - baseline HY distributions

4. `04_datasets_for_models.ipynb` formats datasets to fit models.
5. `05_fit_progression_models.R` fits the models. Runs the same procedure for both LR and CPH.

    **Cohort inclusion.** A cohort is analyzed only if it has at least 20 participants and at least 20 events and 20 non-events. Age and disease duration are centered on the grand mean across all cohorts.

    **Per-cohort estimation.** LR uses Firth-penalized logistic regression, which limits small-sample and separation bias. CPH uses Cox regression with Efron ties.

    **Steps:**

    1. **Base model.** Fits `age + sex + disease duration` in each cohort and combines the results with REML random-effects meta-analysis (`metafor::rma`). These unadjusted associations are the main results.
    2. **HY1-at-baseline model (CPH only).** Fits the base model plus an indicator for HY = 1 at baseline, in cohorts that have this variable. LR skips this step because its event is defined at the same visit.
    3. **Stepwise term selection.** Each candidate term is added to the current model, fit in each cohort and meta-analyzed. A term is kept if the pooled p-value is below 0.05. Terms for selection:
    1. disease duration²
    2. age²
    3. age × sex
    4. disease duration × sex
    5. disease duration³ (considered only if disease duration² was kept)
    6. age³ (considered only if age² was kept)
    4. **Final model: Meta vs. Pooled.** Each term in the selected model is estimated two ways:
    - **Meta:** fit in each cohort, then combine with random-effects meta-analysis.
    - **Pooled:** one fit on all participants together, with cohort as a covariate.

    The Pooled estimate is reported only if it falls inside the Meta 95% CI *and* both methods reach the same conclusion (both significant in the same direction, or both non-significant). Otherwise the Meta estimate is reported. Heterogeneity statistics (τ², I², Cochran's Q) are reported for every term.


6. `06_plot_forest_comparison.R` draws the forest plots.

    - **Main figure:** age, sex and disease duration from the **base model**.
    - **Supplementary figures:** one per stepwise candidate term. A term kept in the final model shows its final-model Meta estimate. A term that was not kept shows its screening-step estimate in grey.
    - **HY1-at-baseline figure:** CPH only, drawn when the HY1-baseline results exist.

---
- Key findings:
    - The impact of medication status on the HY3 outcome is minimal thus making it a stable outcome for modeling
    - Final models:
        - LR: HY3 ~ age + sex + disease duration + age-squared + disease duration-squared (can be pooled)
        - CPH: HY3 ~ age + sex + disease duration (analyze per cohort, then meta-analyze)

---
## Data Statement
Data used in the preparation of this article were obtained from the Global Parkinson’s Genetics Program (GP2; https://gp2.org). All GP2 data are hosted in collaboration with the Accelerating Medicines Partnership in Parkinson’s disease, and are available via application on the website (https://amp-pd.org/register-for-amp-pd). For up-to-date information on GP2 data acquisition, access, and policies, visit https://gp2.org/. Tier 1 data can be accessed by completing a form on the Accelerating Medicines Partnership in Parkinson’s Disease (AMP®-PD) website (https://amp-pd.org/register-for-amp-pd). Tier 2 data access requires approval and a Data Use Agreement signed by your institution. Specifically, we used Tier 2 data from GP2 release 12 (DOI 10.5281/zenodo.20932193). Information on how to access data from the Parkinson's Precision Medicine Initiate (PPMI) can be found by visiting https://www.ppmi-info.org/. 

### Helpful Links

- [GP2 Website](https://gp2.org/)
  - [GP2 Cohort Dashboard](https://gp2.org/cohort-dashboard-advanced/)
- [Introduction to GP2](https://movementdisorders.onlinelibrary.wiley.com/doi/10.1002/mds.28494)
  - [Other GP2 Manuscripts (PubMed)](https://pubmed.ncbi.nlm.nih.gov/?term=%22global+parkinson%27s+genetics+program%22)

## Repository Orientation
- The `analysis/` directory includes all analyses discussed in the manuscript.

<pre> THIS_REPO/ 
  ├── analysis/ 
  |     ├── 01_ppmi_create_datasets.py
  |     ├── 02_hy_ppmi_analysis.R
  |     ├── 03_make_dataset.py
  |     ├── 04_datasets_for_models.ipynb
  |     ├── 05_fit_progression_models.R
  |     └── 06_plot_forest_comparison.R
  ├── LICENSE
  └── README.md 
</pre>

--- 
## Software
Versions listed are those installed in the `python312` conda environment, which provides both Python and R.

| **Software** | **Version(s)** | **Resource URL** | **RRID** | **Notes** |
|--------------|----------------|------------------|----------|-----------|
|Python Programming Language|3.12|http://www.python.org/|RRID:SCR_008394|pandas; numpy; seaborn; matplotlib; scipy; Used for general data wrangling/plotting/analyses|
|R Project for Statistical Computing|v.4.5.3|http://www.r-project.org/|RRID:SCR_001905|tidyverse; survival; metafor; yaml; Used for general data wrangling/plotting/analyses	|

**Python 3.12** (`03_make_dataset.py`)

| Package     | Version | Used for                   |
|-------------|---------|----------------------------|
| pandas      | 2.3.3   | Data loading and wrangling |
| numpy       | 2.4.3   | Numerical operations       |
| matplotlib  | 3.10.9  | Figures                    |
| seaborn     | 0.13.2  | Plot styling               |
| scipy       | 1.17.1  | Smoothing utilities        |

**R 4.5.3** (`05_fit_progression_models.R`, `06_plot_forest_comparison.R`)

| Package   | Version | Used for                                   |
|-----------|---------|--------------------------------------------|
| tidyverse | 2.0.0   | Data wrangling                             |
| survival  | 3.8.6   | Cox proportional hazards models (`coxph`)  |
| metafor   | 5.0.1   | Random-effects meta-analysis (`rma`, REML) |
| yaml      | 2.3.12  | Reading `config.yaml`                      |

The models also use base R: `stats` for logistic regression and the Firth penalized likelihood (fit with `optim`), and `grDevices` for writing the PNG plots.

---
## Citation
If you use this repository or find it helpful for your research, please cite the corresponding manuscript.
