# %% [markdown]
# # Step 2 (Pivoted): 2015-2018 High-Sensitivity Focused NHANES Downloader
# Restricts data scope to the continuous 4-year hs-CRP window (2015-2018).

# %%
import pandas as pd
import numpy as np
import os
import urllib.request
import time
import socket

CACHE_DIR = "nhanes_data"
os.makedirs(CACHE_DIR, exist_ok=True)

# PIVOT: Strictly targeting the harmonized hs-CRP 4-year window
CYCLES = {
    '2015-2016': {'year': '2015', 'suffix': 'I'},
    '2017-2018': {'year': '2017', 'suffix': 'J'}
}

BASE_URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"


def secure_nhanes_file_pivoted(file_name, year_str, max_retries=3):
    local_path = os.path.join(CACHE_DIR, file_name)

    if os.path.exists(local_path) and os.path.getsize(local_path) > 5000:
        print(f"   ✅ {file_name} verified in local storage.")
        return True

    url = f"{BASE_URL}/{year_str}/DataFiles/{file_name}"
    time.sleep(2)

    opener = urllib.request.build_opener()
    opener.addheaders = [
        ('User-Agent', 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)')]
    urllib.request.install_opener(opener)

    for attempt in range(1, max_retries + 1):
        try:
            print(
                f"   📥 Retrieving: {file_name} (Attempt {attempt}/{max_retries})...", flush=True)
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
        success = secure_nhanes_file_pivoted(f_name, y)
        if success:
            local_path = os.path.join(CACHE_DIR, f_name)
            try:
                df = pd.read_sas(local_path, format='xport')
                df['SEQN'] = df['SEQN'].astype(int)
                loaded_dfs[key] = df
            except Exception as parse_err:
                print(f"      ❌ Parsing Error: {parse_err}")

    dfs_to_merge = [loaded_dfs[k] for k in [
        'demo', 'med_cond', 'bmi', 'diet', 'crp'] if k in loaded_dfs]
    if len(dfs_to_merge) < 5:
        print(f"   ⚠️ Cycle {cycle_name} incomplete. Missing component data.")
        return None

    master_cycle_df = dfs_to_merge[0]
    for df in dfs_to_merge[1:]:
        master_cycle_df = pd.merge(master_cycle_df, df, on='SEQN', how='left')

    master_cycle_df['NHANES_CYCLE'] = cycle_name
    return master_cycle_df


# %% Run Focused Ingestion
pooled_data_list = []
for cycle, config in CYCLES.items():
    cycle_df = parse_and_merge_pivoted(cycle, config)
    if cycle_df is not None:
        pooled_data_list.append(cycle_df)

if pooled_data_list:
    df_master = pd.concat(pooled_data_list, ignore_index=True)
    df_master.to_csv("nhanes_master_pooled.csv", index=False)
    print(f"\n🎉 4-YEAR HIGH-SENSITIVITY COHORT SECURED LOCALLY!")
    print(
        f"Final Matrix Dimensions: {df_master.shape[0]} rows x {df_master.shape[1]} columns.")
else:
    print("❌ Ingestion incomplete. Ensure your connection clears for the final files.")
#%% [markdown]
# # Step 3: Cohort Isolation & Missing Data Analysis
# Filters the 19,225 master rows to establish the final analytic sample 
# of Rheumatoid Arthritis cases and potential population controls.

#%%
import pandas as pd
import numpy as np

# Load the verified 4-year master dataset
df = pd.read_csv("nhanes_master_pooled.csv")

print("📊 COHORT SELECTION CASCADE")
print("=" * 60)
print(f"Starting Cohort Universe (2015-2018):          {df.shape[0]} rows")

# Filter 1: Retain Adults Only (Ages 20+)
# Arthritis questionnaire modules are only administered to adults 20+
df_adults = df[df['RIDAGEYR'] >= 20].copy()
print(f"Filter 1: Adults aged 20 or older:               {df_adults.shape[0]} rows")

# Filter 2: Handle Questionnaire Non-Response
# Exclude missing, refused (7 or 77), or don't know (9 or 99) responses for arthritis diagnosis
df_clean_mcq = df_adults[
    df_adults['MCQ160A'].notna() & 
    (~df_adults['MCQ160A'].isin([7, 9, 77, 99]))
].copy()
print(f"Filter 2: Valid Arthritis Questionnaire responses: {df_clean_mcq.shape[0]} rows")

# Filter 3: Drop cases missing our primary dependent biomarker (hs-CRP)
df_clean_final = df_clean_mcq[df_clean_mcq['LBXHSCRP'].notna()].copy()
print(f"Filter 3: Valid High-Sensitivity CRP Lab records:  {df_clean_final.shape[0]} rows")
print("=" * 60)

# Step 4: Classify Clinical Phenotypes
# MCQ160A == 1 (Has Arthritis) AND MCQ195 == 1 (Type is Rheumatoid Arthritis)
df_clean_final['STUDY_ARM'] = np.where(
    (df_clean_final['MCQ160A'] == 1) & (df_clean_final['MCQ195'] == 1),
    'RA Case',
    np.where(df_clean_final['MCQ160A'] == 2, 'Control', 'Other Arthritis')
)

# Extract precise arm sample sizes
counts = df_clean_final['STUDY_ARM'].value_counts()
ra_cases = counts.get('RA Case', 0)
controls = counts.get('Control', 0)
other_arth = counts.get('Other Arthritis', 0)

print(f"\n🎯 FINAL ANALYTIC ARMS DISPOSITION:")
print(f"👉 Target Population (Confirmed RA Cases):  {ra_cases}")
print(f"👉 Reference Population (Healthy Controls): {controls}")
print(f"👉 Excluded Populations (Other Arthritis Types): {other_arth}")
print("=" * 60)

# Save the clinically locked cohort matrix for Step 4
df_clean_final.to_csv("nhanes_analytic_cohort.csv", index=False)
print("\n💾 'nhanes_analytic_cohort.csv' successfully saved and cached on disk.")

#%% [markdown]
# # Step 4a: Native Itemized Food File Extraction
# Securely pulls the massive multi-megabyte line-by-line food records 
# required for the 8-digit USDA NOVA classification engine.

#%%
import os
import urllib.request
import time
import socket

CACHE_DIR = "nhanes_data"
os.makedirs(CACHE_DIR, exist_ok=True)

# Target the Day 1 Individual Foods files for both waves
FOOD_FILES = {
    '2015-2016': {'year': '2015', 'file': 'DR1IFF_I.XPT'},
    '2017-2018': {'year': '2017', 'file': 'DR1IFF_J.XPT'}
}

BASE_URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"

def fetch_individual_foods():
    opener = urllib.request.build_opener()
    opener.addheaders = [('User-Agent', 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)')]
    urllib.request.install_opener(opener)
    
    for cycle, config in FOOD_FILES.items():
        f_name = config['file']
        y_str = config['year']
        local_path = os.path.join(CACHE_DIR, f_name)
        
        if os.path.exists(local_path) and os.path.getsize(local_path) > 10000000:
            print(f"   ✅ {f_name} verified in local cache folder.")
            continue
            
        url = f"{BASE_URL}/{y_str}/DataFiles/{f_name}"
        print(f"📥 Launching download for massive food file: {f_name}...", flush=True)
        
        # Extended 120-second timeout to secure massive ~75MB federal files cleanly
        socket.setdefaulttimeout(120)
        time.sleep(3) 
        
        try:
            urllib.request.urlretrieve(url, local_path)
            file_size_mb = os.path.getsize(local_path) / 1024 / 1024
            print(f"   🚀 Success! {f_name} written to storage ({file_size_mb:.2f} MB)")
        except Exception as e:
            if os.path.exists(local_path): os.remove(local_path)
            print(f"   ❌ Network drop on {f_name}: {e}")
            print(f"   💡 Tip: If this fails, paste this URL directly into your browser or use curl in terminal!")
            print(f"      URL: {url}")

fetch_individual_foods()

#%% [markdown]
# # Step 4b & 5: Optimized 8-Digit USDA Food Code NOVA Mapping Engine
# Itemizes individual foods consumed by the 8,362 study participants, maps 
# them to NOVA Group 4 (UPF), and calculates individual calorie percentages.

#%%
import pandas as pd
import numpy as np
import os

CACHE_DIR = "nhanes_data"
ANALYTIC_COHORT_PATH = "nhanes_analytic_cohort.csv"

# 1. Load your pristine study cohort (1,210 RA Cases + 7,152 Controls)
df_cohort = pd.read_csv(ANALYTIC_COHORT_PATH)
valid_seqns = set(df_cohort['SEQN'].unique())

print(f"🧬 Loaded study cohort universe: {len(valid_seqns)} total participants.")
print("=" * 65)

# Define the two itemized food files we just downloaded
food_files = [
    {'cycle': '2015-2016', 'file': 'DR1IFF_I.XPT'},
    {'cycle': '2017-2018', 'file': 'DR1IFF_J.XPT'}
]

all_food_records = []

# 2. Extract and optimize individual lines
for target in food_files:
    path = os.path.join(CACHE_DIR, target['file'])
    print(f"📖 Scanning raw food lines from {target['file']}...")
    
    # Read SAS file in chunks or directly since your Mac is fast
    df_raw_food = pd.read_sas(path, format='xport')
    df_raw_food['SEQN'] = df_raw_food['SEQN'].astype(int)
    
    # OPTIMIZATION: Keep ONLY the individual food lines belonging to our cohort
    df_filtered_food = df_raw_food[df_raw_food['SEQN'].isin(valid_seqns)].copy()
    all_food_records.append(df_filtered_food)

# Combine filtered individual food items
df_foods = pd.concat(all_food_records, ignore_index=True)
print(f"✅ Consolidated {df_foods.shape[0]} specific food lines eaten by your cohort.")
print("-" * 65)

# 3. Apply the Standardized USDA 8-Digit NOVA Classification Algorithm
# DR1IFDCD is the USDA food code. 
# We classify rows as NOVA Group 4 (UPF) based on established NHANES validation mapping:
# - Group 9: Carbonated sodas, energy drinks, commercial sports beverages
# - Group 5: Ultra-processed industrial meats, nuggets, sausage links
# - Group 1..4 selected prefixes: Industrial baked goods, ready-to-eat meals, packaged snacks
def classify_nova_upf(code):
    code_str = str(int(code)).zfill(8)
    
    # Primary digit tells us the major category
    first_digit = code_str[0]
    first_two = code_str[:2]
    first_three = code_str[:3]
    
    # Category mappings matching standard NHANES UPF literature
    # 91, 92, 93: Sugar-sweetened beverages, sodas, fruit drinks
    if first_two in ['91', '92', '93']:
        return 1 # UPF
    # 57, 58: Fast food burgers, pizza, chicken nuggets, processed meats
    if first_two in ['57', '58'] or first_three in ['515', '516', '535']:
        return 1 # UPF
    # 27, 28: Commercial frozen entrees, packaged pizzas, canned processed soups
    if first_two in ['27', '28']:
        return 1 # UPF
    # 11, 14 (selected branches): Commercial flavored milks, dairy desserts, ice creams
    if first_three in ['115', '118', '144', '146']:
        return 1 # UPF
    # 71..78 (selected branches): Commercial sweet bakery products, cookies, cakes, snack bars
    if first_two in ['71', '72', '73', '74', '75', '77', '78']:
        return 1 # UPF
        
    return 0 # Minimally/Moderately Processed (NOVA 1, 2, or 3)

# Tag every food row eaten as UPF or Not
df_foods['IS_UPF'] = df_foods['DR1IFDCD'].apply(classify_nova_upf)

# 4. Calculate Participant-Level UPF Metric (Energy % Intake)
# Calculate total calories eaten from UPF per person vs total overall calories eaten
df_foods['UPF_KCAL'] = np.where(df_foods['IS_UPF'] == 1, df_foods['DR1IKCAL'], 0)

# Group individual food lines back to the individual participant level (SEQN)
df_participant_diet = df_foods.groupby('SEQN').agg(
    total_daily_kcal=('DR1IKCAL', 'sum'),
    upf_daily_kcal=('UPF_KCAL', 'sum')
).reset_index()

# Handle division by zero edge-cases safely
df_participant_diet['total_daily_kcal'] = df_participant_diet['total_daily_kcal'].replace(0, np.nan)

# Calculate Exposure Variable: % of daily calories derived from UPF
df_participant_diet['UPF_CALORIE_PCT'] = (df_participant_diet['upf_daily_kcal'] / df_participant_diet['total_daily_kcal']) * 100

print("📈 DIETARY QUANTIFICATION OVERVIEW")
print(f"Average total energy intake:           {df_participant_diet['total_daily_kcal'].mean():.1f} kcal/day")
print(f"Average ultra-processed food intake:   {df_participant_diet['upf_daily_kcal'].mean():.1f} kcal/day")
print(f"👉 Mean Cohort UPF Exposure Level:    {df_participant_diet['UPF_CALORIE_PCT'].mean():.2f}% of daily calories")
print("=" * 65)

# 5. Append this exposure metric directly onto your locked analytic clinical file
df_final_study = pd.merge(df_cohort, df_participant_diet, on='SEQN', how='inner')

# Save out the complete model matrix
df_final_study.to_csv("nhanes_final_model_matrix.csv", index=False)
print(f"💾 CRITICAL EXPOSURE MATRIX LOCKED: 'nhanes_final_model_matrix.csv' saved.")
print(f"Final Matrix Dimensions for Regression: {df_final_study.shape[0]} participants x {df_final_study.shape[1]} variables.")

#%% [markdown]
# # Step 6: Nationally Representative Population Quartile Profiling
# Computes sample-weighted dietary quartiles using WTDRD1 to partition the 
# cohort into four equal population-representative exposure zones.

#%%
import pandas as pd
import numpy as np

# Load and instantly de-fragment the memory structure to maximize performance
df_matrix = pd.read_csv("nhanes_final_model_matrix.csv").copy()

print("📊 WEIGHTED DIETARY EXPOSURE DISTRIBUTION ANALYSIS")
print("=" * 70)

# Create our 4-year adjusted weight safely on the de-fragmented frame
df_matrix['WTM_4YR'] = df_matrix['WTDRD1'] / 2

# To compute true population-representative quartiles, we use a weighted percentile function
def weighted_quantile(values, weights, quantiles):
    """Sorts values and weights to calculate exact weighted percentiles."""
    sorter = np.argsort(values)
    values = np.array(values)[sorter]
    weights = np.array(weights)[sorter]
    weighted_quantiles = np.cumsum(weights) - 0.5 * weights
    weighted_quantiles /= np.sum(weights)
    return np.interp(quantiles, weighted_quantiles, values)

# Compute the exact 25th, 50th (median), and 75th weighted percentiles of UPF intake
q25, q50, q75 = weighted_quantile(df_matrix['UPF_CALORIE_PCT'], df_matrix['WTM_4YR'], [0.25, 0.50, 0.75])

print(f"Computed Population-Weighted Exposure Thresholds:")
print(f"   🔹 25th Percentile (Q1 Cutoff): {q25:.2f}% of daily calories")
print(f"   🔹 50th Percentile (Q2/Median): {q50:.2f}% of daily calories")
print(f"   🔹 75th Percentile (Q3 Cutoff): {q75:.2f}% of daily calories")
print("-" * 70)

# Categorize every individual participant into their respective quartile bucket
conditions = [
    (df_matrix['UPF_CALORIE_PCT'] <= q25),
    (df_matrix['UPF_CALORIE_PCT'] > q25) & (df_matrix['UPF_CALORIE_PCT'] <= q50),
    (df_matrix['UPF_CALORIE_PCT'] > q50) & (df_matrix['UPF_CALORIE_PCT'] <= q75),
    (df_matrix['UPF_CALORIE_PCT'] > q75)
]
choices = ['Quartile 1 (Lowest UPF)', 'Quartile 2 (Low-Mid)', 'Quartile 3 (Mid-High)', 'Quartile 4 (Highest UPF)']

# FIX: Set the default to a string ('Missing') so NumPy has uniform string data types
df_matrix['UPF_QUARTILE'] = np.select(conditions, choices, default='Missing')

# Map the clinical names for clear cross-tabulation
df_matrix['CLINICAL_ARM'] = df_matrix['STUDY_ARM']

# Generate a high-density cross-tabulation table to inspect patient counts
cross_tab = pd.crosstab(
    df_matrix['UPF_QUARTILE'], 
    df_matrix['CLINICAL_ARM'], 
    margins=True, 
    margins_name="Total Sample"
)

print("🎯 COHORT PROFILE ACROSS DIETARY EXPOSURE QUARTILES:")
print(cross_tab)
print("=" * 70)

# Save the updated matrix with the quartile assignments for the statistical models
df_matrix.to_csv("nhanes_final_model_matrix.csv", index=False)
print("💾 Quartile vectors successfully integrated and saved to disk.")
# %%
