# Ultra-processed food intake and inflammation in rheumatoid arthritis

Cross-sectional analysis of NHANES 2015–2018 data examining whether higher ultra-processed food (UPF) intake is associated with systemic inflammation (high-sensitivity CRP) among adults with self-reported rheumatoid arthritis (RA).

**Primary question:** Among adults with RA, is greater UPF intake (as a share of daily calories) associated with higher CRP?

## Repository layout

| Path | Description |
|------|-------------|
| `main.py` | End-to-end pipeline: download, merge, cohort definition, NOVA mapping, quartiles, survey-weighted GEE models, Table 1 output |
| `plan.md` | Step-by-step implementation checklist |
| `nhanes_data/` | Cached NHANES SAS transport (`.XPT`) files |
| `nhanes_master_pooled.csv` | Pooled 2015–2018 demographic, clinical, dietary, and lab data |
| `nhanes_analytic_cohort.csv` | Adults with valid arthritis and hs-CRP records; RA / control / other arms |
| `nhanes_final_model_matrix.csv` | Analysis-ready file with UPF exposure, quartiles, survey design vars |

## Requirements

- Python 3.10+
- [Git LFS](https://git-lfs.com/) (two large dietary recall files are stored with LFS)

## Setup

```bash
git clone git@github.com:ahint01/upf-crp.git
cd upf-crp
git lfs pull

python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Running the analysis

`main.py` is a linear script with `# %%` cell markers (VS Code / Cursor “Run Cell” or Jupyter-compatible runners). Run the full file:

```bash
python main.py
```

On first run, smaller NHANES files are downloaded from the CDC into `nhanes_data/`. Large individual dietary recall files (`DR1IFF_I.XPT`, `DR1IFF_J.XPT`) are included via Git LFS. Intermediate CSVs are written at each stage and reused on subsequent runs.

## Pipeline overview

1. **Ingest** — Download and merge DEMO, MCQ, BMX, DR1TOT, and HSCRP for 2015–2016 and 2017–2018.
2. **Cohort** — Adults ≥20 with valid arthritis questionnaire and hs-CRP; classify RA cases (`MCQ160A` + `MCQ195`), controls, and other arthritis.
3. **Diet / NOVA** — Line-level dietary recalls; rule-based 8-digit USDA code mapping to UPF; participant % calories from UPF.
4. **Quartiles** — Survey-weighted (`WTDRD1`) UPF quartiles within the analytic sample.
5. **Descriptives** — CRP distributions by UPF quartile in the RA cohort.
6. **Survey design** — Attach `SDMVSTRA` and `SDMVPSU` from DEMO files.
7. **Regression** — GEE on log(CRP) with nested cluster IDs (`SDMVSTRA` + `SDMVPSU`), survey weights, and covariates (age, sex, race/ethnicity, PIR, BMI, total kcal). RA subset with hs-CRP ≤ 10 mg/L.
8. **Table 1** — Descriptive characteristics of the RA cohort by UPF quartile.

See `plan.md` for the detailed checklist.

## Key variables

| Role | Variable |
|------|----------|
| Exposure | % daily calories from NOVA Group 4 (UPF), quartiled |
| Outcome | hs-CRP (`LBXHSCRP`), log-transformed in models |
| Population | Self-reported RA (`STUDY_ARM == 'RA Case'`) |
| Covariates | Age, sex, race/ethnicity, BMI, poverty–income ratio, total energy intake |
| Weights | `WTDRD1` / 2 for 4-year pooled dietary weights |

## Data sources

Public-use NHANES files from the [CDC NHANES site](https://wwwn.cdc.gov/nchs/nhanes/):

- **2015–2016** (suffix `I`): 2015 release year in download URLs  
- **2017–2018** (suffix `J`): 2017 release year in download URLs  

This repository is for reproducible research. NHANES data are de-identified public-use files; analyses use survey weights and complex design variables where noted. Results are exploratory and not clinical guidance.

## License

Add a license file if you plan to share or publish this code. NHANES data remain subject to [NCHS data use policies](https://www.cdc.gov/nchs/nhanes/index.htm).
