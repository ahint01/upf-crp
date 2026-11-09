# Ultra-processed food intake and inflammation in rheumatoid arthritis

Analysis code for the poster **"Ultra-Processed Food Intake and Inflammatory Burden Among Adults with Rheumatoid Arthritis: A NHANES Analysis"**.

Julia Tripodo, B.S., Andrew Hinton, B.S., Naziha Slimani, M.D.
College of Medicine, University of Central Florida

Cross-sectional analysis of pooled NHANES 2015–2018 data asking whether, among adults with self-reported rheumatoid arthritis (RA), a higher share of calories from ultra-processed foods (UPF) is associated with higher high-sensitivity C-reactive protein (hs-CRP).

## Results

| | |
|---|---|
| Analytic sample | 9,199 adults ≥20 with a valid arthritis response, hs-CRP and a day-1 dietary recall |
| RA participants | 1,126 |
| Adjusted model (RA, hs-CRP ≤10 mg/L, complete covariates) | 886 |
| Mean share of calories from UPF (whole sample) | 47.4% |
| Median hs-CRP in RA, lowest vs highest UPF quartile | 2.20 vs 2.81 mg/L |

Adjusted hs-CRP ratios versus the lowest UPF quartile (survey-weighted GEE, adjusted for age, sex, race/ethnicity, poverty-income ratio, BMI and total energy intake):

| UPF quartile | hs-CRP ratio | 95% CI | p |
|---|---|---|---|
| Q2 (low-mid) | 1.26 | 1.00–1.58 | 0.055 |
| Q3 (mid-high) | 1.28 | 1.04–1.59 | 0.022 |
| Q4 (highest) | 1.30 | 0.99–1.70 | 0.058 |

These are cross-sectional associations, not causal estimates. RA status is self-reported and diet comes from a single 24-hour recall. Results are exploratory and not clinical guidance.

## Running the analysis

Requires Python 3.10+.

```bash
git clone git@github.com:ahint01/upf-crp.git
cd upf-crp

python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt

python main.py
```

`main.py` is a linear script with `# %%` cell markers, so it also runs cell by cell in VS Code, Cursor or Jupyter.

The first run downloads the NHANES public-use files from the CDC into `nhanes_data/` (the dietary recall files are large, so allow a few minutes) and writes intermediate CSVs to the repo root. Later runs reuse the cached files. Data files are not tracked in git; everything is regenerated from the CDC source.

## Pipeline

1. **Ingest:** download and merge DEMO, MCQ, BMX, DR1TOT and HSCRP for 2015–2016 (`_I`) and 2017–2018 (`_J`), the two cycles with hs-CRP.
2. **Cohort:** adults ≥20 with a valid arthritis answer and hs-CRP. RA = ever told arthritis (`MCQ160A == 1`) and type rheumatoid (`MCQ195 == 1`).
3. **UPF exposure:** classify each food in the day-1 individual foods file (`DR1IFF`) by its 8-digit USDA food code using an adapted NOVA approach, then compute each participant's percent of calories from UPF.
4. **Quartiles:** weighted UPF quartiles (`WTDRD1 / 2`, the NCHS rule for pooling two cycles) on the full analytic sample.
5. **Descriptives:** hs-CRP median, IQR and share >3 mg/L by quartile in RA.
6. **Survey design:** attach strata (`SDMVSTRA`) and PSUs (`SDMVPSU`).
7. **Model:** GEE on log(hs-CRP) with sampling weights and stratum–PSU clusters, in RA participants with hs-CRP ≤10 mg/L to exclude likely acute inflammation.
8. **Table 1:** RA characteristics by UPF quartile.

## Data

Public-use files from the [CDC NHANES site](https://wwwn.cdc.gov/nchs/nhanes/), subject to [NCHS data use policies](https://www.cdc.gov/nchs/nhanes/index.htm).
