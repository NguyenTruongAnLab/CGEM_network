#!/usr/bin/env python3
"""
Compare C-GEM model output with March 2025 validation data
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Setup paths
case_dir = Path('INPUT/Cases/Mekong_Delta_Full')
output_dir = Path('OUTPUT/Mekong_Delta_Full/CSV')
fig_dir = Path('OUTPUT/Mekong_Delta_Full/figures')
fig_dir.mkdir(parents=True, exist_ok=True)

# Load validation data
val_df = pd.read_csv(case_dir / 'Validation' / 'Mekong_4_branches_Mar_2025.csv')

# Convert units to model units
# O2: mg/L -> umol/L (x31.25)
# NO3, NH4: mgN/L -> umolN/L (x71.43)
# DOC: mgC/L -> umolC/L (x83.33)
# CH4: ugC/L -> nmol/L (x83.33)
# N2O: ugN/L -> nmol/L (x35.71)
val_df['O2_umol'] = val_df['O2 (mg/L)'] * 31.25
val_df['NO3_umol'] = val_df['NO3-N (mgN/L)'] * 71.43
val_df['NH4_umol'] = val_df['NH4 (mgN/L)'] * 71.43
val_df['DOC_umol'] = val_df['DOC (mgC/L)'] * 83.33
val_df['CH4_nmol'] = val_df['ugC-CH4/L'] / 12 * 1000
val_df['N2O_nmol'] = val_df['ugN-N2O/L'] / 28 * 1000
val_df['Alk_umol'] = val_df['Alkalinity (mmol/L)'] * 1000

# Branch mapping
branch_map = {
    'Hau_River': 'Hau_River',
    'Ham_Luong': 'Ham_Luong', 
    'Co_Chien': 'Co_Chien',
    'My_Tho': 'My_Tho'
}

def load_model_output(branch, variable):
    """Load model CSV output for a branch and variable"""
    fname = output_dir / f"{branch}_{variable}.csv"
    if not fname.exists():
        return None, None
    
    df = pd.read_csv(fname)
    # Get last row (final timestep)
    last_row = df.iloc[-1]
    
    # Parse distance columns
    distances = []
    values = []
    for col in df.columns:
        if col == 'Time_s':
            continue
        if 'km' in col:
            try:
                dist = float(col.replace('km', ''))
                distances.append(dist)
                values.append(last_row[col])
            except:
                pass
    
    return np.array(distances), np.array(values)

# Create comparison figures
fig, axes = plt.subplots(4, 4, figsize=(20, 16))
plt.suptitle('C-GEM Model vs March 2025 Validation Data', fontsize=14, y=1.02)

variables = [
    ('salinity', 'Surface Salinity', 'Salinity [PSU]', 1.0),
    ('o2', 'O2_umol', 'O2 [µmol/L]', 1.0),
    ('spm', 'TSS (mg/L)', 'SPM [mg/L]', 1.0),
    ('pco2', 'pCO2 (ppm)', 'pCO2 [ppm]', 1.0),
]

for row, branch in enumerate(['Hau_River', 'Ham_Luong', 'Co_Chien', 'My_Tho']):
    # Get validation data for this branch
    bval = val_df[val_df['River'] == branch].sort_values('Distance')
    
    for col, (var_model, var_val, ylabel, scale) in enumerate(variables):
        ax = axes[row, col]
        
        # Load model output
        dist_m, val_m = load_model_output(branch, var_model)
        
        if dist_m is not None:
            ax.plot(dist_m, val_m * scale, 'b-', linewidth=2, label='Model')
        
        # Plot validation data
        if var_val in bval.columns:
            ax.scatter(bval['Distance'], bval[var_val], c='red', s=50, 
                      marker='o', label='Obs Mar 2025', zorder=5)
        
        ax.set_xlabel('Distance from mouth [km]')
        ax.set_ylabel(ylabel)
        ax.set_title(f'{branch}')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=8)
        
        # Set axis limits based on variable
        if var_model == 'salinity':
            ax.set_ylim(0, 35)
        elif var_model == 'o2':
            ax.set_ylim(100, 300)
        elif var_model == 'spm':
            ax.set_ylim(0, 50)
        elif var_model == 'pco2':
            ax.set_ylim(0, 6000)

plt.tight_layout()
plt.savefig(fig_dir / 'model_vs_validation.png', dpi=150, bbox_inches='tight')
print(f"Saved: {fig_dir / 'model_vs_validation.png'}")

# Create second figure for nutrients and GHGs
fig2, axes2 = plt.subplots(4, 4, figsize=(20, 16))
plt.suptitle('C-GEM Model vs March 2025 - Biogeochemistry', fontsize=14, y=1.02)

variables2 = [
    ('no3', 'NO3_umol', 'NO3 [µmol/L]', 1.0),
    ('ph', 'pH', 'pH', 1.0),
    ('ch4', 'CH4_nmol', 'CH4 [nmol/L]', 1.0),
    ('n2o', 'N2O_nmol', 'N2O [nmol/L]', 1.0),
]

for row, branch in enumerate(['Hau_River', 'Ham_Luong', 'Co_Chien', 'My_Tho']):
    bval = val_df[val_df['River'] == branch].sort_values('Distance')
    
    for col, (var_model, var_val, ylabel, scale) in enumerate(variables2):
        ax = axes2[row, col]
        
        dist_m, val_m = load_model_output(branch, var_model)
        
        if dist_m is not None:
            ax.plot(dist_m, val_m * scale, 'b-', linewidth=2, label='Model')
        
        if var_val in bval.columns:
            ax.scatter(bval['Distance'], bval[var_val], c='red', s=50,
                      marker='o', label='Obs Mar 2025', zorder=5)
        
        ax.set_xlabel('Distance from mouth [km]')
        ax.set_ylabel(ylabel)
        ax.set_title(f'{branch}')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=8)
        
        if var_model == 'no3':
            ax.set_ylim(0, 80)
        elif var_model == 'ph':
            ax.set_ylim(7.0, 8.5)
        elif var_model == 'ch4':
            ax.set_ylim(0, 300)
        elif var_model == 'n2o':
            ax.set_ylim(0, 200)

plt.tight_layout()
plt.savefig(fig_dir / 'model_vs_validation_biogeo.png', dpi=150, bbox_inches='tight')
print(f"Saved: {fig_dir / 'model_vs_validation_biogeo.png'}")

# Calculate RMSE for each variable
print("\n=== MODEL VALIDATION STATISTICS ===\n")
print(f"{'Branch':<15} {'Variable':<12} {'RMSE':<10} {'R²':<10} {'Bias':<10}")
print("-" * 60)

for branch in ['Hau_River', 'Ham_Luong', 'Co_Chien', 'My_Tho']:
    bval = val_df[val_df['River'] == branch].sort_values('Distance')
    
    for var_model, var_val, ylabel, scale in variables + variables2:
        dist_m, val_m = load_model_output(branch, var_model)
        
        if dist_m is None or var_val not in bval.columns:
            continue
        
        # Interpolate model to observation distances
        obs_dist = bval['Distance'].values
        obs_val = bval[var_val].dropna()
        
        if len(obs_val) < 3:
            continue
        
        # Simple interpolation
        model_interp = np.interp(obs_dist[:len(obs_val)], dist_m, val_m * scale)
        
        # Calculate stats
        diff = model_interp - obs_val.values
        rmse = np.sqrt(np.mean(diff**2))
        bias = np.mean(diff)
        
        # R² (if enough points)
        if len(obs_val) > 2:
            ss_res = np.sum(diff**2)
            ss_tot = np.sum((obs_val.values - np.mean(obs_val.values))**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        else:
            r2 = np.nan
        
        print(f"{branch:<15} {var_model:<12} {rmse:<10.2f} {r2:<10.2f} {bias:<+10.2f}")

print("\n" + "=" * 60)
print("Note: Negative bias = model underestimates, positive = overestimates")
