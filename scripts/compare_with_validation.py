#!/usr/bin/env python3
"""
Compare C-GEM model output with March 2025 validation data
Generates comparison figures and calculates error metrics
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Paths
case_dir = Path('INPUT/Cases/Mekong_Delta_Full')
output_dir = Path('OUTPUT/Mekong_Delta_Full/CSV')
fig_dir = Path('OUTPUT/Mekong_Delta_Full/figures')
fig_dir.mkdir(parents=True, exist_ok=True)

# Load validation data (in model units)
val_file = case_dir / 'Validation' / 'Mekong_Mar2025_model_units.csv'
val_df = pd.read_csv(val_file)
print(f"Loaded {len(val_df)} validation points from {val_file}")

# Branch mapping
branches = ['Hau_River', 'Ham_Luong', 'Co_Chien', 'My_Tho']
branch_colors = {
    'Hau_River': '#1f77b4',
    'Ham_Luong': '#ff7f0e', 
    'Co_Chien': '#2ca02c',
    'My_Tho': '#d62728'
}

def load_model_output(branch, variable):
    """Load model output for a branch and variable"""
    filepath = output_dir / f"{branch}_{variable}.csv"
    if not filepath.exists():
        return None, None
    
    df = pd.read_csv(filepath)
    # Get last timestep (final state)
    last_row = df.iloc[-1]
    
    # Extract distances and values
    distances = []
    values = []
    for col in df.columns[1:]:  # Skip Time_s
        if 'km' in col:
            dist = float(col.replace('km', ''))
            distances.append(dist)
            values.append(last_row[col])
    
    # Apply unit conversions
    # Model outputs CH4/N2O in µmol/L (µM), validation is in nmol/L
    # 1 µM = 1000 nmol/L
    if variable in ['ch4', 'n2o']:
        values = [v * 1000 for v in values]  # Convert µM → nmol/L
    
    return np.array(distances), np.array(values)

def calculate_metrics(obs, mod):
    """Calculate error metrics"""
    if len(obs) == 0 or len(mod) == 0:
        return {'n': 0}
    
    # Only use non-NaN pairs
    mask = ~np.isnan(obs) & ~np.isnan(mod)
    obs = obs[mask]
    mod = mod[mask]
    
    if len(obs) < 3:
        return {'n': len(obs)}
    
    bias = np.mean(mod - obs)
    rmse = np.sqrt(np.mean((mod - obs)**2))
    mae = np.mean(np.abs(mod - obs))
    
    # R² (coefficient of determination)
    if np.std(obs) > 0:
        r, p = stats.pearsonr(obs, mod)
        r2 = r**2
    else:
        r2 = np.nan
    
    return {
        'n': len(obs),
        'bias': bias,
        'rmse': rmse,
        'mae': mae,
        'r2': r2,
        'obs_mean': np.mean(obs),
        'mod_mean': np.mean(mod)
    }

# Variables to compare
variables = [
    ('salinity', 'Salinity_PSU', 'Salinity [PSU]', 0, 35),
    ('o2', 'O2_umol_L', 'O2 [µmol/L]', 100, 300),
    ('spm', 'SPM_mg_L', 'SPM [mg/L]', 0, 50),
    ('pco2', 'pCO2_ppm', 'pCO2 [ppm]', 0, 6000),
    ('ph', 'pH', 'pH', 7.0, 8.5),
    ('no3', 'NO3_umol_L', 'NO3 [µmol/L]', 0, 80),
    ('nh4', 'NH4_umol_L', 'NH4 [µmol/L]', 0, 10),
    ('toc', 'TOC_umol_L', 'TOC [µmol/L]', 0, 250),
    ('at', 'Alkalinity_umol_L', 'Alkalinity [µmol/L]', 1000, 2500),
    ('ch4', 'CH4_nmol_L', 'CH4 [nmol/L]', 0, 300),
    ('n2o', 'N2O_nmol_L', 'N2O [nmol/L]', 0, 200),
]

# Create main comparison figure
fig, axes = plt.subplots(4, 3, figsize=(18, 20))
axes = axes.flatten()

all_metrics = {}

for idx, (model_var, val_var, ylabel, ymin, ymax) in enumerate(variables):
    ax = axes[idx]
    
    all_metrics[model_var] = {}
    
    for branch in branches:
        # Get validation data for this branch
        bval = val_df[val_df['Branch'] == branch].sort_values('Distance_km')
        if len(bval) == 0:
            continue
        
        # Get model output
        model_dist, model_vals = load_model_output(branch, model_var)
        if model_dist is None:
            continue
        
        # Plot model output
        ax.plot(model_dist, model_vals, '-', color=branch_colors[branch], 
                linewidth=2, alpha=0.8, label=f'{branch} (model)')
        
        # Plot validation data
        ax.scatter(bval['Distance_km'], bval[val_var], 
                   color=branch_colors[branch], s=60, edgecolors='black',
                   linewidths=1, zorder=5, marker='o')
        
        # Interpolate model to validation points for metrics
        if len(model_dist) > 0 and len(bval) > 0:
            model_interp = np.interp(bval['Distance_km'], model_dist, model_vals)
            metrics = calculate_metrics(bval[val_var].values, model_interp)
            all_metrics[model_var][branch] = metrics
    
    ax.set_xlabel('Distance from mouth (km)')
    ax.set_ylabel(ylabel)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0, 110)
    ax.grid(True, alpha=0.3)
    ax.set_title(ylabel)
    
    if idx == 0:
        ax.legend(loc='upper right', fontsize=8)

# Remove empty subplot
axes[-1].axis('off')

plt.suptitle('C-GEM Model vs March 2025 Field Data', fontsize=16, y=1.01)
plt.tight_layout()
plt.savefig(fig_dir / 'model_validation_comparison.png', dpi=300, bbox_inches='tight')
print(f"Saved: {fig_dir / 'model_validation_comparison.png'}")

# Print metrics summary
print("\n" + "=" * 80)
print("MODEL VALIDATION METRICS")
print("=" * 80)

for var in all_metrics:
    print(f"\n{var.upper()}:")
    for branch in all_metrics[var]:
        m = all_metrics[var][branch]
        if m['n'] >= 3:
            print(f"  {branch:15s}: n={m['n']:2d}, RMSE={m['rmse']:7.2f}, Bias={m['bias']:+7.2f}, R²={m['r2']:.3f}")

# Create scatter plots for key variables
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

key_vars = [('salinity', 'Salinity_PSU', 'Salinity [PSU]'),
            ('o2', 'O2_umol_L', 'O2 [µmol/L]'),
            ('pco2', 'pCO2_ppm', 'pCO2 [ppm]'),
            ('ph', 'pH', 'pH'),
            ('spm', 'SPM_mg_L', 'SPM [mg/L]'),
            ('toc', 'TOC_umol_L', 'TOC [µmol/L]')]

for idx, (model_var, val_var, ylabel) in enumerate(key_vars):
    ax = axes[idx]
    
    all_obs = []
    all_mod = []
    
    for branch in branches:
        bval = val_df[val_df['Branch'] == branch].sort_values('Distance_km')
        if len(bval) == 0:
            continue
        
        model_dist, model_vals = load_model_output(branch, model_var)
        if model_dist is None:
            continue
        
        # Interpolate
        model_interp = np.interp(bval['Distance_km'], model_dist, model_vals)
        
        ax.scatter(bval[val_var], model_interp, 
                   color=branch_colors[branch], s=40, alpha=0.7, 
                   label=branch, edgecolors='white', linewidths=0.5)
        
        all_obs.extend(bval[val_var].values)
        all_mod.extend(model_interp)
    
    # 1:1 line
    all_obs = np.array(all_obs)
    all_mod = np.array(all_mod)
    mask = ~np.isnan(all_obs) & ~np.isnan(all_mod)
    
    if np.sum(mask) > 0:
        vmin = min(np.nanmin(all_obs[mask]), np.nanmin(all_mod[mask]))
        vmax = max(np.nanmax(all_obs[mask]), np.nanmax(all_mod[mask]))
        ax.plot([vmin, vmax], [vmin, vmax], 'k--', linewidth=1, label='1:1')
        
        # Calculate overall metrics
        metrics = calculate_metrics(all_obs[mask], all_mod[mask])
        if metrics['n'] >= 3:
            ax.text(0.05, 0.95, f"RMSE={metrics['rmse']:.1f}\nR²={metrics['r2']:.2f}", 
                    transform=ax.transAxes, fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_xlabel(f'Observed {ylabel}')
    ax.set_ylabel(f'Modeled {ylabel}')
    ax.set_title(ylabel)
    ax.grid(True, alpha=0.3)
    
    if idx == 0:
        ax.legend(loc='lower right', fontsize=8)

plt.suptitle('Model vs Observation Scatter Plots', fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(fig_dir / 'model_validation_scatter.png', dpi=300, bbox_inches='tight')
print(f"Saved: {fig_dir / 'model_validation_scatter.png'}")

# Create detailed profile comparison for each branch
for branch in branches:
    fig, axes = plt.subplots(3, 4, figsize=(20, 12))
    axes = axes.flatten()
    
    bval = val_df[val_df['Branch'] == branch].sort_values('Distance_km')
    
    plot_vars = [
        ('salinity', 'Salinity_PSU', 'Salinity [PSU]'),
        ('o2', 'O2_umol_L', 'O2 [µmol/L]'),
        ('pco2', 'pCO2_ppm', 'pCO2 [ppm]'),
        ('ph', 'pH', 'pH'),
        ('spm', 'SPM_mg_L', 'SPM [mg/L]'),
        ('no3', 'NO3_umol_L', 'NO3 [µmol/L]'),
        ('nh4', 'NH4_umol_L', 'NH4 [µmol/L]'),
        ('po4', 'PO4_umol_L', 'PO4 [µmol/L]'),
        ('toc', 'TOC_umol_L', 'TOC [µmol/L]'),
        ('at', 'Alkalinity_umol_L', 'Alk [µmol/L]'),
        ('ch4', 'CH4_nmol_L', 'CH4 [nmol/L]'),
        ('n2o', 'N2O_nmol_L', 'N2O [nmol/L]'),
    ]
    
    for idx, (model_var, val_var, ylabel) in enumerate(plot_vars):
        ax = axes[idx]
        
        model_dist, model_vals = load_model_output(branch, model_var)
        
        if model_dist is not None:
            ax.plot(model_dist, model_vals, '-', color='blue', 
                    linewidth=2, label='Model')
        
        if val_var in bval.columns:
            ax.scatter(bval['Distance_km'], bval[val_var], 
                       color='red', s=80, edgecolors='black',
                       linewidths=1, zorder=5, label='Observed')
        
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=8)
    
    plt.suptitle(f'{branch} - Model vs Validation', fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(fig_dir / f'validation_{branch}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {fig_dir / 'validation_{branch}.png'}")

plt.close('all')

print("\n" + "=" * 80)
print("VALIDATION COMPLETE")
print("=" * 80)
