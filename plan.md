# Project Plan: Ultra-Processed Food Intake and Inflammatory Burden in Rheumatoid Arthritis (NHANES)

## 📌 Project Overview
- **Goal:** Evaluate whether higher ultra-processed food (UPF) intake is associated with increased systemic inflammation (CRP) among adults with self-reported rheumatoid arthritis (RA).
- **Primary Question:** Among adults with RA, is greater UPF intake associated with higher CRP levels?
- **Design:** Cross-sectional analysis using NHANES data.

## 📊 Key Variables
* **Exposure:** UPF intake (% of total daily calories from NOVA Group 4 foods).
* **Outcome:** CRP level (continuous log-transformed and/or dichotomized as elevated CRP).
* **Population:** Adults ages 20+ years with self-reported RA.
* **Covariates:** Age, sex, race/ethnicity, BMI, smoking status, diabetes status, physical activity, total caloric intake, and socioeconomic status (PIR).

---

## 🗺️ Step-by-Step Implementation Tracker

- [x] **Step 1:** Identify NHANES cycles including dietary recall, CRP, RA status, and covariates.
- [x] **Step 2:** Download and merge NHANES files (Demographics, Dietary, Labs, Questionnaires) using Python.
- [x] **Step 3:** Define the RA cohort based on self-reported arthritis questionnaire variables.
- [x] **Step 4:** Map NHANES 8-digit food codes to the NOVA classification system (UPF vs. non-UPF).
- [x] **Step 5:** Calculate participant-level percentage of daily calories derived from UPFs.
- [x] **Step 6:** Categorize UPF intake into quartiles and verify continuous distributions.
- [/] **Step 7:** Perform descriptive statistics and compare baseline CRP across UPF groups.
- [ ] **Step 8:** Run multivariable survey-weighted regression models adjusting for covariates.
- [ ] **Step 9:** Conduct sensitivity analyses (e.g., excluding CRP > 10 mg/L to rule out acute infection).
- [ ] **Step 10:** Format and export results/tables for ACR Convergence abstract submission.