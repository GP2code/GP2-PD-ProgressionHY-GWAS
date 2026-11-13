import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import UnivariateSpline
from scipy.stats import gaussian_kde
import os
import sys
import logging
from datetime import datetime

# ── Setup logging ─────────────────────────────────────────────────────────────
os.makedirs('logs', exist_ok=True)

log_file = 'logs/make_dataset_log.txt'
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s',
    handlers=[
        logging.FileHandler(log_file, mode='w'),
        logging.StreamHandler(sys.stdout)
    ]
)

logging.info(f"Dataset Preparation Log - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
logging.info("=" * 80)

# ── Load ──────────────────────────────────────────────────────────────────────
## Data from GP2 R12 release
logging.info("\nLoading GP2 R12 HY data (cross-sectional and longitudinal cohorts)...")
dcs = pd.read_csv('r12_hy/hy_cs_GP2_2026-08-10.csv') # cross-sectional data for checking duplicates and cohort design
dlt = pd.read_csv('r12_hy/hy_lt_GP2_2026-08-10.csv') # longitudinal data

dt = pd.concat([dlt, dcs]).drop_duplicates(subset=['GP2ID', 'visit_month'], keep='first')
dt['cohort'] = dt.study
dt['IID'] = dt['GP2ID']
dt['y'] = dt['hoehn_and_yahr_stage']
dt['SEX'] = dt['biological_sex_for_qc'].map({'Male': 1, 'Female': 2})
dt['disease_duration_new'] = dt['disease_duration']
dt['study_days'] = dt['visit_month'] 
# remove the study n<20
logging.info("\nRemoving cohorts with <20 patients...")
cohort_counts = dt.groupby('cohort')['IID'].nunique()
cohorts_to_remove = cohort_counts[cohort_counts < 20].index.tolist()
dt = dt[~dt['cohort'].isin(cohorts_to_remove)]
logging.info(f"Remaining cohorts: {dt['cohort'].unique()}")
logging.info(f"\nLoaded: {dt.shape[0]} rows, {dt['cohort'].nunique()} cohorts")
logging.info(f"Unique individuals: {dt['IID'].nunique()}")

# ── 1. Disease duration at current visit (years) ──────────────────────────────
dt['disease_duration'] = dt['disease_duration_new'] + (dt['study_days'] / 12)
logging.info("\n[1] Calculating disease duration at current visit...")
logging.info(f"   Disease duration range: {dt['disease_duration'].min():.2f} - {dt['disease_duration'].max():.2f} years")
logging.info(f"   Mean: {dt['disease_duration'].mean():.2f} ± {dt['disease_duration'].std():.2f} years")

# ── 2. Age at onset ───────────────────────────────────────────────────────────
dt['age_at_onset'] = dt['age_at_baseline'] - dt['disease_duration_new']
logging.info("\n[2] Calculating age at onset...")
logging.info(f"   Age at onset range: {dt['age_at_onset'].min():.1f} - {dt['age_at_onset'].max():.1f} years")
logging.info(f"   Mean: {dt['age_at_onset'].mean():.1f} ± {dt['age_at_onset'].std():.1f} years")

# ── 3. Age at onset (mean-centered) ──────────────────────────────────────────
dt['age_at_onset_c'] = dt['age_at_onset'] - dt['age_at_onset'].mean()
logging.info(f"\n[3] Mean-centering age at onset (mean = {dt['age_at_onset'].mean():.1f} years)")

# ── 4. Sex label ──────────────────────────────────────────────────────────────
dt['sex'] = dt['SEX'].map({1: 'Male', 2: 'Female'})
logging.info("\n[4] sex distribution:")
logging.info(dt['sex'].value_counts().to_string())

# ── 5. Rename columns ─────────────────────────────────────────────────────────
dt = dt.rename(columns={'IID': 'id', 'y': 'HY'})
logging.info(f"\n[5] columns after rename: {dt.columns.tolist()}")

# ── 6. Left-truncation flag (lt) per cohort ───────────────────────────────────
# Longitudinal cohort  → any ID appears ≥2 times → lt = 1
# Cross-sectional cohort → all IDs appear exactly once → lt = 0
cohort_lt = (
    dt.groupby('cohort')['id']
    .apply(lambda x: int((x.value_counts() >= 2).any()))
    .rename('lt')
)
dt = dt.join(cohort_lt, on='cohort')

logging.info("\n[6] Cohort design assessment:")
summary = (
    dt.groupby('cohort')
    .agg(n_rows=('id', 'count'), n_unique_ids=('id', 'nunique'), lt=('lt', 'first'))
    .assign(design=lambda x: x['lt'].map({1: 'longitudinal', 0: 'cross-sectional'}))
)
logging.info(summary.to_string())

# ── 7. Select final columns ───────────────────────────────────────────────────
final_cols = ['id', 'HY', 'disease_duration', 'sex', 'age_at_onset', 'age_at_onset_c', 
              'age_at_baseline', 'cohort', 'lt']
if 'ancestry' in dt.columns:
    final_cols.append('ancestry')

df_all = dt[final_cols].copy()
logging.info(f"\n[7] Final dataset (all data): {df_all.shape[0]} rows x {df_all.shape[1]} cols")
logging.info(f"   Unique individuals: {df_all['id'].nunique()}")
logging.info(f"   Cohorts: {df_all['cohort'].nunique()}")

# ── 8. Data Cleaning for Analysis ─────────────────────────────────────────────
# Apply cleaning filters: 
# - Remove disease_duration > 20 years
# - Remove HY == 0 (abnormal values)
# - Collapse intermediate HY stages (1.5→1, 2.5→2, 3.5→3, 4.5→4)
logging.info("\n[8] Applying data cleaning filters...")
logging.info(f"   Before cleaning: {df_all.shape[0]} rows")

df_analysis = df_all.copy()
df_analysis = df_analysis[df_analysis['disease_duration'] <= 20]
logging.info(f"   After disease_duration ≤ 20: {df_analysis.shape[0]} rows") 

df_analysis = df_analysis[df_analysis['HY'] != 0]
logging.info(f"   After removing HY=0: {df_analysis.shape[0]} rows")

# Remove HY4+ within 1 years of disease duration (likely end-stage patients with atypical parkinsonism)
df_analysis = df_analysis[~((df_analysis['HY'] >= 4) & (df_analysis['disease_duration'] <= 1))]
logging.info(f"   Remaining rows after removing HY4+ within 1 year: {df_analysis.shape[0]}")

# Collapse intermediate HY stages
df_analysis['HY'] = df_analysis['HY'].replace({1.5: 1, 2.5: 2, 3.5: 3, 4.5: 4})
logging.info(f"   HY distribution after collapsing intermediate stages:")
logging.info(df_analysis['HY'].value_counts().sort_index().to_string())

# Final sex counts
logging.info("\n[4] sex distribution:")
logging.info(df_analysis['sex'].value_counts().to_string())

# ── 9. Save both datasets ─────────────────────────────────────────────────────
os.makedirs('temp', exist_ok=True)
df_all.to_csv('temp/hy_analysis_dataset_all.csv', index=False)
logging.info(f"\n[9] Saved: temp/hy_analysis_dataset_all.csv")
logging.info(f"     Observations: {df_all.shape[0]}, Individuals: {df_all['id'].nunique()}, Cohorts: {df_all['cohort'].nunique()}")

df_analysis.to_csv('temp/hy_analysis_dataset.csv', index=False)
logging.info(f"\n     Saved: temp/hy_analysis_dataset.csv")
logging.info(f"     Observations: {df_analysis.shape[0]}, Individuals: {df_analysis['id'].nunique()}, Cohorts: {df_analysis['cohort'].nunique()}")
logging.info(f"     Removed: {df_all.shape[0] - df_analysis.shape[0]} rows ({100*(df_all.shape[0] - df_analysis.shape[0])/df_all.shape[0]:.1f}%)")

# ══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION & REPORTING
# ══════════════════════════════════════════════════════════════════════════════

logging.info("\n" + "=" * 80)
logging.info("Generating figures and tables for both datasets...")

# Set matplotlib style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# Create separate report directories for each dataset
os.makedirs('report/data_summary_all', exist_ok=True)
os.makedirs('report/data_summary_clean', exist_ok=True)

# Loop through both datasets
datasets = [
    {'name': 'all', 'df': df_all, 'report_dir': 'report/data_summary_all'},
    {'name': 'clean', 'df': df_analysis, 'report_dir': 'report/data_summary_clean'}
]

for dataset_info in datasets:
    dataset_name = dataset_info['name']
    dt = dataset_info['df']
    report_dir = dataset_info['report_dir']
    
    logging.info("\n" + "-" * 80)
    logging.info(f"Processing dataset: {dataset_name} ({dt.shape[0]} rows)")
    logging.info("-" * 80)
    
    # Prepare baseline data (first observation per patient)
    dt_baseline = dt.sort_values(['id', 'disease_duration']).groupby('id').first().reset_index()
    
    # ── Figure A: Disease Duration Histogram per HY Score ────────────────────────
    logging.info(f"\n[{dataset_name}] Figure A: Creating disease duration histogram by HY score...")
    
    hy_values = sorted(dt['HY'].dropna().unique())
    n_hy = len(hy_values)
    ncols = min(3, n_hy)
    nrows = int(np.ceil(n_hy / ncols))
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    if n_hy == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for idx, hy in enumerate(hy_values):
        ax = axes[idx]
        data = dt[dt['HY'] == hy]['disease_duration'].dropna()
        ax.hist(data, bins=30, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Disease Duration (years)', fontsize=10)
        ax.set_ylabel('Frequency', fontsize=10)
        ax.set_title(f'HY = {hy} (n={len(data)})', fontsize=11, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
    
    # Hide unused subplots
    for idx in range(n_hy, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    fig_a_path = f'{report_dir}/fig_disease_duration_by_hy.png'
    plt.savefig(fig_a_path, dpi=300, bbox_inches='tight')
    plt.close()
    logging.info(f"Saved: {fig_a_path}")
    
    # ── Figure B: HY Score vs Disease Duration with Spline (by Sex) ──────────────
    logging.info(f"\n[{dataset_name}] Figure B: Creating HY by disease duration with spline (stratified by sex)...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    for idx, sex_val in enumerate(['Male', 'Female']):
        ax = axes[idx]
        data_sex = dt[dt['sex'] == sex_val].dropna(subset=['HY', 'disease_duration'])
        
        # Scatter plot with jitter to show overlapping discrete HY values
        jitter = np.random.normal(0, 0.08, size=len(data_sex))  # Small vertical jitter
        ax.scatter(data_sex['disease_duration'], data_sex['HY'] + jitter, 
                   alpha=0.3, s=10, label='Observations')
        
        # Rolling window mean with CI
        if len(data_sex) > 50:
            # Sort data by disease duration
            sorted_data = data_sex.sort_values('disease_duration')
            x_data = sorted_data['disease_duration'].values
            y_data = sorted_data['HY'].values
            
            try:
                # Define rolling window parameters
                window_size = max(50, int(len(data_sex) * 0.1))  # Use 10% of data as window
                
                # Calculate rolling mean
                x_smooth = []
                y_smooth = []
                y_lower = []
                y_upper = []
                
                for i in range(len(x_data)):
                    # Define window around current point
                    half_window = window_size // 2
                    start_idx = max(0, i - half_window)
                    end_idx = min(len(x_data), i + half_window + 1)
                    
                    window_y = y_data[start_idx:end_idx]
                    
                    # Calculate mean and confidence interval
                    mean_y = np.mean(window_y)
                    std_y = np.std(window_y)
                    sem_y = std_y / np.sqrt(len(window_y))
                    
                    x_smooth.append(x_data[i])
                    y_smooth.append(mean_y)
                    y_lower.append(mean_y - 1.96 * sem_y)  # 95% CI
                    y_upper.append(mean_y + 1.96 * sem_y)
                
                # Subsample for smoother plot (every 5th point)
                step = max(1, len(x_smooth) // 200)
                x_smooth = np.array(x_smooth[::step])
                y_smooth = np.clip(np.array(y_smooth[::step]), 0, 5)  # Clip to valid HY range
                y_lower = np.clip(np.array(y_lower[::step]), 0, 5)
                y_upper = np.clip(np.array(y_upper[::step]), 0, 5)
                
                ax.plot(x_smooth, y_smooth, 'r-', linewidth=2, label='Rolling mean')
                ax.fill_between(x_smooth, y_lower, y_upper, alpha=0.2, color='red', 
                                label='95% CI')
            except Exception as e:
                logging.warning(f"Rolling mean smoothing failed for {sex_val}: {e}")
        
        ax.set_xlabel('Disease Duration (years)', fontsize=11)
        ax.set_ylabel('HY Score', fontsize=11)
        ax.set_title(f'{sex_val} (n={len(data_sex)})', fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=9)
        ax.grid(alpha=0.3)
    
    plt.tight_layout()
    fig_b_path = f'{report_dir}/fig_hy_by_duration_sex.png'
    plt.savefig(fig_b_path, dpi=300, bbox_inches='tight')
    plt.close()
    logging.info(f"Saved: {fig_b_path}")
    
    # ── Figure D: Age at Onset Distribution by Cohort ────────────────────────────
    logging.info(f"\n[{dataset_name}] Figure D: Creating age at onset distribution by cohort...")
    
    cohorts = sorted(dt_baseline['cohort'].unique())
    n_cohorts = len(cohorts)
    
    fig, ax = plt.subplots(figsize=(max(12, n_cohorts*0.8), 6))
    data_list = [dt_baseline[dt_baseline['cohort'] == c]['age_at_onset'].dropna().values 
                 for c in cohorts]
    bp = ax.boxplot(data_list, labels=cohorts, patch_artist=True)
    
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
        patch.set_alpha(0.7)
    
    ax.set_xlabel('Cohort', fontsize=11)
    ax.set_ylabel('Age at Onset (years)', fontsize=11)
    ax.set_title('Age at Onset Distribution by Cohort', fontsize=13, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    fig_d_path = f'{report_dir}/fig_age_at_onset_by_cohort.png'
    plt.savefig(fig_d_path, dpi=300, bbox_inches='tight')
    plt.close()
    logging.info(f"Saved: {fig_d_path}")
    
    # ── Figure E: Sample Sizes per Cohort ────────────────────────────────────────
    logging.info(f"\n[{dataset_name}] Figure E: Creating sample sizes per cohort...")
    
    cohort_stats = dt.groupby('cohort').agg(
        n_patients=('id', 'nunique'),
        n_observations=('id', 'count')
    ).reset_index()
    
    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(cohort_stats))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, cohort_stats['n_patients'], width, 
                   label='Unique Patients', alpha=0.8, color='steelblue')
    bars2 = ax.bar(x + width/2, cohort_stats['n_observations'], width, 
                   label='Total Observations', alpha=0.8, color='coral')
    
    ax.set_xlabel('Cohort', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title('Sample Sizes per Cohort', fontsize=13, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(cohort_stats['cohort'], rotation=45, ha='right')
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    
    fig_e_path = f'{report_dir}/fig_sample_sizes.png'
    plt.savefig(fig_e_path, dpi=300, bbox_inches='tight')
    plt.close()
    logging.info(f"Saved: {fig_e_path}")
    
    # ── Figure F: HY Distribution by Cohort (Grid) ───────────────────────────────
    logging.info(f"\n[{dataset_name}] Figure F: Creating HY distribution by cohort (grid layout)...")
    
    n_cohorts = len(cohorts)
    ncols_f = min(4, n_cohorts)
    nrows_f = int(np.ceil(n_cohorts / ncols_f))
    
    fig, axes = plt.subplots(nrows_f, ncols_f, figsize=(4*ncols_f, 3.5*nrows_f))
    if n_cohorts == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for idx, cohort in enumerate(cohorts):
        ax = axes[idx]
        data_cohort = dt_baseline[dt_baseline['cohort'] == cohort]
        hy_counts = data_cohort['HY'].value_counts().sort_index()
        
        hy_vals = hy_counts.index.tolist()
        counts = hy_counts.values.tolist()
        
        ax.bar(range(len(hy_vals)), counts, alpha=0.7, edgecolor='black')
        ax.set_xticks(range(len(hy_vals)))
        ax.set_xticklabels([f'{v}' for v in hy_vals], fontsize=9)
        ax.set_xlabel('HY Score', fontsize=9)
        ax.set_ylabel('Count (Baseline)', fontsize=9)
        ax.set_title(f'{cohort} (n={len(data_cohort)})', fontsize=10, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
    
    # Hide unused subplots
    for idx in range(n_cohorts, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    fig_f_path = f'{report_dir}/fig_hy_distribution_by_cohort.png'
    plt.savefig(fig_f_path, dpi=300, bbox_inches='tight')
    plt.close()
    logging.info(f"Saved: {fig_f_path}")
    
    # ── Table 1: Cohort Baseline Summary ─────────────────────────────────────────
    logging.info(f"\n[{dataset_name}] Table 1: Creating cohort baseline summary table...")
    
    cohorts_sorted = sorted(dt_baseline['cohort'].unique())
    table_data = {}
    
    for cohort in cohorts_sorted:
        data_cohort = dt_baseline[dt_baseline['cohort'] == cohort]
        data_cohort_full = dt[dt['cohort'] == cohort]
        
        # Basic counts
        n_patients = len(data_cohort)
        n_obs = len(data_cohort_full)
        
        # Sex distribution
        n_male = (data_cohort['sex'] == 'Male').sum()
        n_female = (data_cohort['sex'] == 'Female').sum()
        pct_male = 100 * n_male / n_patients if n_patients > 0 else 0
        pct_female = 100 * n_female / n_patients if n_patients > 0 else 0
        
        # Age at baseline
        age_baseline_mean = data_cohort['age_at_baseline'].mean()
        age_baseline_sd = data_cohort['age_at_baseline'].std()
        
        # Age at onset (unconverted from centered)
        age_onset_mean = data_cohort['age_at_onset'].mean()
        age_onset_sd = data_cohort['age_at_onset'].std()
        
        # Disease duration at baseline
        dd_mean = data_cohort['disease_duration'].mean()
        dd_sd = data_cohort['disease_duration'].std()
        dd_median = data_cohort['disease_duration'].median()
        dd_q1 = data_cohort['disease_duration'].quantile(0.25)
        dd_q3 = data_cohort['disease_duration'].quantile(0.75)
        
        # Follow-up duration (for longitudinal cohorts only)
        is_longitudinal = data_cohort['lt'].iloc[0] == 1 if len(data_cohort) > 0 else False
        if is_longitudinal:
            followup_per_patient = data_cohort_full.groupby('id')['disease_duration'].apply(
                lambda x: x.max() - x.min() if len(x) > 1 else 0
            )
            fu_mean = followup_per_patient.mean()
            fu_sd = followup_per_patient.std()
            fu_median = followup_per_patient.median()
            fu_q1 = followup_per_patient.quantile(0.25)
            fu_q3 = followup_per_patient.quantile(0.75)
        else:
            fu_mean = fu_sd = fu_median = fu_q1 = fu_q3 = np.nan
        
        # HY distribution at baseline
        hy_counts = data_cohort['HY'].value_counts()
        hy_data = {}
        for hy_val in sorted(dt_baseline['HY'].dropna().unique()):
            count = hy_counts.get(hy_val, 0)
            pct = 100 * count / n_patients if n_patients > 0 else 0
            hy_data[f'HY={hy_val}'] = f'{count} ({pct:.1f}%)'
        
        # Compile cohort column
        table_data[cohort] = {
            'n_patients': n_patients,
            'n_observations': n_obs,
            'Male': f'{n_male} ({pct_male:.1f}%)',
            'Female': f'{n_female} ({pct_female:.1f}%)',
            'Age at baseline (mean±SD)': f'{age_baseline_mean:.1f} ± {age_baseline_sd:.1f}',
            'Age at onset (mean±SD)': f'{age_onset_mean:.1f} ± {age_onset_sd:.1f}',
            'Disease duration baseline (mean±SD)': f'{dd_mean:.2f} ± {dd_sd:.2f}',
            'Disease duration baseline (median[Q1,Q3])': f'{dd_median:.2f} [{dd_q1:.2f}, {dd_q3:.2f}]',
            'Follow-up duration (mean±SD)': f'{fu_mean:.2f} ± {fu_sd:.2f}' if not np.isnan(fu_mean) else 'N/A (cross-sectional)',
            'Follow-up duration (median[Q1,Q3])': f'{fu_median:.2f} [{fu_q1:.2f}, {fu_q3:.2f}]' if not np.isnan(fu_median) else 'N/A (cross-sectional)',
            **hy_data
        }
    
    # Create DataFrame
    table1 = pd.DataFrame(table_data)
    table1_path = f'{report_dir}/table1_cohort_baseline_summary.csv'
    table1.to_csv(table1_path)
    logging.info(f"Saved: {table1_path}")
    logging.info(f"\nTable 1 Preview ({dataset_name}):")
    logging.info(table1.to_string())

logging.info("\n" + "=" * 80)
logging.info("All figures and tables generated successfully for both datasets!")
logging.info(f"Log file saved: {log_file}")
logging.info("=" * 80)

