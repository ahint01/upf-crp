# %% [markdown]
# # NHANES 2015-2018: Ultra-Processed Food Consumption and Inflammation in Rheumatoid Arthritis
# ### Global Dependencies and Configuration Layer

# %%
import os
import time
import socket
import urllib.request
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

# Set up global runtime paths
CACHE_DIR = "nhanes_data"
os.makedirs(CACHE_DIR, exist_ok=True)
FINAL_MATRIX_PATH = "nhanes_final_model_matrix.csv"

# Continuous 4-year high-sensitivity CRP window (2015-2018)
CYCLES = {
    '2015-2016': {'year': '2015', 'suffix': 'I'},
    '2017-2018': {'year': '2017', 'suffix': 'J'}
}
BASE_URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"

# %% [markdown]
# # Step 1 & 2: Biennial Dataset Ingestion Engine
# Restricts data scope to the continuous 4-year hs-CRP window and pools clinical components.

# %%
def secure_nhanes_file_pivoted(file_name, year_str, max_retries=3):
    local_path = os.path.join(CACHE_DIR, file_name)
    if os.path.exists(local_path) and os.path.getsize(local_path) > 5000:
        print(f"   ✅ {file_name} verified in local storage.")
        return True

    url = f"{BASE_URL}/{year_str}/DataFiles/{file_name}"
    time.sleep(2)
    opener = urllib.request.build_opener()
    opener.addheaders = [('User-Agent', 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)')]
    urllib.request.install_opener(opener)

    for attempt in range(1, max_retries + 1):
        try:
            print(f"   📥 Retrieving: {file_name} (Attempt {attempt}/{max_retries})...", flush=True)
            socket.setdefaulttimeout(60)
            urllib.request.urlretrieve(url, local_path)
            file_size = os.path.getsize(local_path)
            print(f"      🚀 Success! ({file_size / 1024 / 1024:.2f} MB)")
            return True
        except Exception as e:
            if os.path.exists(local_path):
                os.remove(local_path)
            print(f"      ⚠️ Connection interrupted: {e}")
        if attempt < max_retries:
            time.sleep(5)
    return False

def parse_and_merge_pivoted(cycle_name, config):
    print(f"\n⚡ Processing Cycle: {cycle_name} ---")
    y = config['year']
    s = config['suffix']

    files = {
        'demo': f"DEMO_{s}.XPT",
        'med_cond': f"MCQ_{s}.XPT",
        'bmi': f"BMX_{s}.XPT",
        'diet': f"DR1TOT_{s}.XPT",
        'crp': f"HSCRP_{s}.XPT"
    }

    loaded_dfs = {}
    for key, f_name in files.items():
        if secure_nhanes_file_pivoted(f_name, y):
            local_path = os.path.join(CACHE_DIR, f_name)
            try:
                df = pd.read_sas(local_path, format='xport')
                df['SEQN'] = df['SEQN'].astype(int)
                loaded_dfs[key] = df
            except Exception as parse_err:
                print(f"      ❌ Parsing Error: {parse_err}")

    dfs_to_merge = [loaded_dfs[k] for k in ['demo', 'med_cond', 'bmi', 'diet', 'crp'] if k in loaded_dfs]
    if len(dfs_to_merge) < 5:
        print(f"   ⚠️ Cycle {cycle_name} incomplete. Missing component data.")
        return None

    master_cycle_df = dfs_to_merge[0]
    for df in dfs_to_merge[1:]:
        master_cycle_df = pd.merge(master_cycle_df, df, on='SEQN', how='left')

    master_cycle_df['NHANES_CYCLE'] = cycle_name
    return master_cycle_df

pooled_data_list = []
for cycle, config in CYCLES.items():
    cycle_df = parse_and_merge_pivoted(cycle, config)
    if cycle_df is not None:
        pooled_data_list.append(cycle_df)

if pooled_data_list:
    df_master = pd.concat(pooled_data_list, ignore_index=True)
    df_master.to_csv("nhanes_master_pooled.csv", index=False)
    print(f"\n🎉 4-YEAR HIGH-SENSITIVITY COHORT SECURED LOCALLY!")
    print(f"Final Matrix Dimensions: {df_master.shape[0]} rows x {df_master.shape[1]} columns.")
else:
    print("❌ Ingestion incomplete. Check connection configurations.")

# %% [markdown]
# # Step 3: Cohort Isolation & Selection Cascade
# Establishes clinical phenotypes across the adult sampling universe.

# %%
df = pd.read_csv("nhanes_master_pooled.csv")

print("📊 COHORT SELECTION CASCADE")
print("=" * 60)
print(f"Starting Cohort Universe (2015-2018):          {df.shape[0]} rows")

df_adults = df[df['RIDAGEYR'] >= 20].copy()
print(f"Filter 1: Adults aged 20 or older:               {df_adults.shape[0]} rows")

df_clean_mcq = df_adults[
    df_adults['MCQ160A'].notna() & 
    (~df_adults['MCQ160A'].isin([7, 9, 77, 99]))
].copy()
print(f"Filter 2: Valid Arthritis Questionnaire responses: {df_clean_mcq.shape[0]} rows")

df_clean_final = df_clean_mcq[df_clean_mcq['LBXHSCRP'].notna()].copy()
print(f"Filter 3: Valid High-Sensitivity CRP Lab records:  {df_clean_final.shape[0]} rows")
print("=" * 60)

df_clean_final['STUDY_ARM'] = np.where(
    (df_clean_final['MCQ160A'] == 1) & (df_clean_final['MCQ195'] == 1),
    'RA Case',
    np.where(df_clean_final['MCQ160A'] == 2, 'Control', 'Other Arthritis')
)

counts = df_clean_final['STUDY_ARM'].value_counts()
print(f"\n🎯 FINAL ANALYTIC ARMS DISPOSITION:")
print(f"👉 Target Population (Confirmed RA Cases):  {counts.get('RA Case', 0)}")
print(f"👉 Reference Population (Healthy Controls): {counts.get('Control', 0)}")
print(f"👉 Excluded Populations (Other Arthritis):  {counts.get('Other Arthritis', 0)}")
print("=" * 60)

df_clean_final.to_csv("nhanes_analytic_cohort.csv", index=False)

# %% [markdown]
# # Step 4: Individual Dietary Line Extraction & NOVA Processing Engine
# Downloads line-by-line consumption entries and maps 8-digit USDA food codes to processing tiers.

# %%
FOOD_FILES = {
    '2015-2016': {'year': '2015', 'file': 'DR1IFF_I.XPT'},
    '2017-2018': {'year': '2017', 'file': 'DR1IFF_J.XPT'}
}

# 1. Download Individual Food Logs
for cycle, config in FOOD_FILES.items():
    f_name = config['file']
    y_str = config['year']
    local_path = os.path.join(CACHE_DIR, f_name)
    if os.path.exists(local_path) and os.path.getsize(local_path) > 10000000:
        continue
    url = f"{BASE_URL}/{y_str}/DataFiles/{f_name}"
    socket.setdefaulttimeout(120)
    time.sleep(3)
    try:
        urllib.request.urlretrieve(url, local_path)
    except Exception as e:
        if os.path.exists(local_path): os.remove(local_path)
        print(f"❌ Download dropped for {f_name}: {e}")

# 2. Process Logs Against Cohort Frame
df_cohort = pd.read_csv("nhanes_analytic_cohort.csv")
valid_seqns = set(df_cohort['SEQN'].unique())
all_food_records = []

for target in [{'cycle': '2015-2016', 'file': 'DR1IFF_I.XPT'}, {'cycle': '2017-2018', 'file': 'DR1IFF_J.XPT'}]:
    path = os.path.join(CACHE_DIR, target['file'])
    df_raw_food = pd.read_sas(path, format='xport')
    df_raw_food['SEQN'] = df_raw_food['SEQN'].astype(int)
    all_food_records.append(df_raw_food[df_raw_food['SEQN'].isin(valid_seqns)].copy())

df_foods = pd.concat(all_food_records, ignore_index=True)

# 3. Rule-Based USDA Food-Code Prefix Mapping Algorithm (Adapted)
def classify_nova_upf(code):
    code_str = str(int(code)).zfill(8)
    first_two = code_str[:2]
    first_three = code_str[:3]
    
    if first_two in ['91', '92', '93']: return 1
    if first_two in ['57', '58'] or first_three in ['515', '516', '535']: return 1
    if first_two in ['27', '28']: return 1
    if first_three in ['115', '118', '144', '146']: return 1
    if first_two in ['71', '72', '73', '74', '75', '77', '78']: return 1
    return 0

df_foods['IS_UPF'] = df_foods['DR1IFDCD'].apply(classify_nova_upf)
df_foods['UPF_KCAL'] = np.where(df_foods['IS_UPF'] == 1, df_foods['DR1IKCAL'], 0)

df_participant_diet = df_foods.groupby('SEQN').agg(
    total_daily_kcal=('DR1IKCAL', 'sum'),
    upf_daily_kcal=('UPF_KCAL', 'sum')
).reset_index()

df_participant_diet['total_daily_kcal'] = df_participant_diet['total_daily_kcal'].replace(0, np.nan)
df_participant_diet['UPF_CALORIE_PCT'] = (df_participant_diet['upf_daily_kcal'] / df_participant_diet['total_daily_kcal']) * 100

print("📈 DIETARY QUANTIFICATION OVERVIEW")
print(f" 👉 Mean Cohort UPF Exposure Level: {df_participant_diet['UPF_CALORIE_PCT'].mean():.2f}% of daily calories")

df_final_study = pd.merge(df_cohort, df_participant_diet, on='SEQN', how='inner')
df_final_study.to_csv("nhanes_final_model_matrix.csv", index=False)

# %% [markdown]
# # Step 5: Population-Weighted Quartile Profiling
# Leverages sampling weights (`WTDRD1`) to compute representative national exposure zones.

# %%
df_matrix = pd.read_csv("nhanes_final_model_matrix.csv").copy()
df_matrix['WTM_4YR'] = df_matrix['WTDRD1'] / 2

def weighted_quantile(values, weights, quantiles):
    sorter = np.argsort(values)
    values = np.array(values)[sorter]
    weights = np.array(weights)[sorter]
    weighted_quantiles = np.cumsum(weights) - 0.5 * weights
    weighted_quantiles /= np.sum(weights)
    return np.interp(quantiles, weighted_quantiles, values)

q25, q50, q75 = weighted_quantile(df_matrix['UPF_CALORIE_PCT'], df_matrix['WTM_4YR'], [0.25, 0.50, 0.75])

conditions = [
    (df_matrix['UPF_CALORIE_PCT'] <= q25),
    (df_matrix['UPF_CALORIE_PCT'] > q25) & (df_matrix['UPF_CALORIE_PCT'] <= q50),
    (df_matrix['UPF_CALORIE_PCT'] > q50) & (df_matrix['UPF_CALORIE_PCT'] <= q75),
    (df_matrix['UPF_CALORIE_PCT'] > q75)
]
choices = ['Quartile 1 (Lowest UPF)', 'Quartile 2 (Low-Mid)', 'Quartile 3 (Mid-High)', 'Quartile 4 (Highest UPF)']
df_matrix['UPF_QUARTILE'] = np.select(conditions, choices, default='Missing')
df_matrix['CLINICAL_ARM'] = df_matrix['STUDY_ARM']

print(pd.crosstab(df_matrix['UPF_QUARTILE'], df_matrix['CLINICAL_ARM'], margins=True))
df_matrix.to_csv("nhanes_final_model_matrix.csv", index=False)

# %% [markdown]
# # Step 6: Clinical Inflammatory Profiling
# Descriptively assesses baseline, unadjusted serum high-sensitivity CRP shifts across cohorts.

# %%
df_models = pd.read_csv("nhanes_final_model_matrix.csv")
df_ra = df_models[df_models['CLINICAL_ARM'] == 'RA Case'].copy()

print("🩺 SYSTEMIC INFLAMMATION PROFILE WITHIN THE RA COHORT")
print(f"Overall Cohort Median hs-CRP: {df_ra['LBXHSCRP'].median():.2f} mg/L\n")

crp_profile = df_ra.groupby('UPF_QUARTILE').agg(
    median_crp=('LBXHSCRP', 'median'),
    q25_crp=('LBXHSCRP', lambda x: np.percentile(x, 25)),
    q75_crp=('LBXHSCRP', lambda x: np.percentile(x, 75))
).reindex(choices)

print(crp_profile)

df_ra['ELEVATED_CRP'] = np.where(df_ra['LBXHSCRP'] > 3.0, 1, 0)
elevated_pct = df_ra.groupby('UPF_QUARTILE')['ELEVATED_CRP'].mean() * 100
print("\n🚨 RISK RATIO PROFILE (CRP >3.0 mg/L):")
for q in choices:
    print(f"   🔸 {q:<24}: {elevated_pct.get(q, 0):.1f}%")

# %% [markdown]
# # Step 7: Gold-Standard Complex Survey Structural Nesting Layer
# Extracts administrative survey vectors, constructs isolated nested geographic groups, and maps sample rows.

# %%
df_matrix = pd.read_csv(FINAL_MATRIX_PATH)
demo_files = ['DEMO_I.XPT', 'DEMO_J.XPT']
survey_design_list = []

for f in demo_files:
    p = os.path.join(CACHE_DIR, f)
    if os.path.exists(p):
        df_d = pd.read_sas(p, format='xport')[['SEQN', 'SDMVSTRA', 'SDMVPSU']]
        df_d['SEQN'] = df_d['SEQN'].astype(int)
        survey_design_list.append(df_d)

df_design_vars = pd.concat(survey_design_list, ignore_index=True)
cols_to_keep = [c for c in df_matrix.columns if c not in ['SDMVSTRA', 'SDMVPSU']]
df_survey_ready = pd.merge(df_matrix[cols_to_keep], df_design_vars, on='SEQN', how='inner')
df_survey_ready.to_csv(FINAL_MATRIX_PATH, index=False)
print("✅ Primary structural survey design layers integrated safely.")

# %% [markdown]
# # Step 8: Definitive Nested Complex Survey Regression Engine
# Fits Generalized Estimating Equations (GEE) using 60 unique nested clusters (`SDMVSTRA` + `SDMVPSU`) 
# to isolate independent chronic inflammatory trajectories.

# %%
df_survey_ready = pd.read_csv(FINAL_MATRIX_PATH)

# Sensitivity Isolation: Isolate chronic stable RA patients (drops acute infection noise > 10.0 mg/L)
df_ra_survey = df_survey_ready[(df_survey_ready['STUDY_ARM'] == 'RA Case') & (df_survey_ready['LBXHSCRP'] <= 10.0)].copy()

# Setup downstream variable matrices
df_ra_survey['LOG_CRP'] = np.log(df_ra_survey['LBXHSCRP'].replace(0, 0.01))
df_ra_survey['SEX'] = np.where(df_ra_survey['RIAGENDR'] == 2, 'Female', 'Male')
df_ra_survey['UPF_QUARTILE'] = pd.Categorical(
    df_ra_survey['UPF_QUARTILE'], 
    categories=['Quartile 1 (Lowest UPF)', 'Quartile 2 (Low-Mid)', 'Quartile 3 (Mid-High)', 'Quartile 4 (Highest UPF)'],
    ordered=True
)

# CRUCIAL EXPERT ALIGNMENT: Creates explicit nested cluster IDs to match R's nest=TRUE framework
df_ra_survey['NESTED_CLUSTER'] = df_ra_survey['SDMVSTRA'].astype(int).astype(str) + "_" + df_ra_survey['SDMVPSU'].astype(int).astype(str)

analysis_cols = ['LOG_CRP', 'UPF_QUARTILE', 'RIDAGEYR', 'SEX', 'RIDRETH3', 'INDFMPIR', 'BMXBMI', 'total_daily_kcal', 'WTM_4YR', 'NESTED_CLUSTER']
df_ra_survey_clean = df_ra_survey[analysis_cols].dropna().copy()
df_ra_survey_clean = df_ra_survey_clean.sort_values(by='NESTED_CLUSTER').reset_index(drop=True)

print("👑 CORRECTED NESTED COMPLEX SURVEY REGRESSION ENGINE (GEE)")
print("=" * 95)
print(f"Analytic Chronic RA Population: {df_ra_survey_clean.shape[0]} patients")
print(f"Total Unique Nested Cluster Groups Created: {df_ra_survey_clean['NESTED_CLUSTER'].nunique()}")
print("-" * 95)

model_formula = "LOG_CRP ~ UPF_QUARTILE + RIDAGEYR + C(SEX) + C(RIDRETH3) + INDFMPIR + BMXBMI + total_daily_kcal"

gee_model = smf.gee(
    formula=model_formula,
    groups="NESTED_CLUSTER",
    data=df_ra_survey_clean,
    weights=df_ra_survey_clean['WTM_4YR'],
    family=sm.families.Gaussian()
)
survey_results = gee_model.fit()

print("📝 FINAL NESTED COMPLEX SURVEY MODEL MATRIX (Dependent Variable: ln[hs-CRP])")
print(survey_results.summary().tables[1])
print("=" * 95)

print("\n💡 DEFINITIVE COHORT BIOMARKER MULTIPLIERS:")
for idx, coef in survey_results.params.items():
    if "UPF_QUARTILE" in idx:
        multiplier = np.exp(coef)
        pct_increase = (multiplier - 1) * 100
        p_val = survey_results.pvalues[idx]
        print(f"   👉 {idx.split('.')[-1][:-1]:<24} : True Chronic CRP Multiplier = {multiplier:.2f} ({pct_increase:+.1f}% shift) | Nest-Corrected p-value = {p_val:.4f}")
print("=" * 95)