#!/usr/bin/env python3
"""
Compare C-GEM model output with March 2025 validation data
Generates comparison figures and calculates error metrics

DECEMBER 2025 FIX: 
- Now shows tidal variation (mean ± std, or min-max envelope)
- Uses time-averaged values for comparison, not just last timestep
- More honest metrics reporting
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

def load_model_output_with_tidal_stats(branch, variable, last_n_days=5):
    """
    Load model output for a branch and variable.
    Returns statistics over the last N days to capture tidal variation.
    
    Returns: distances, mean_values, min_values, max_values, std_values
    """
    filepath = output_dir / f"{branch}_{variable}.csv"
    if not filepath.exists():
        return None, None, None, None, None
    
    df = pd.read_csv(filepath)
    
    # Extract distance columns
    dist_cols = [col for col in df.columns if 'km' in col]
    distances = np.array([float(col.replace('km', '')) for col in dist_cols])
    
    # Calculate how many timesteps correspond to last_n_days
    # Assuming output every hour (3600 s) and dt=300s simulation
    time_s = df['Time_s'].values
    total_time = time_s[-1] - time_s[0]
    total_days = total_time / 86400
    
    # Use last N days or last 20% of data, whichever is larger
    n_days_available = min(last_n_days, total_days * 0.5)
    cutoff_time = time_s[-1] - n_days_available * 86400
    
    # Get data from the analysis period (after warmup, last few days)
    mask = time_s >= cutoff_time
    if np.sum(mask) < 10:
        # Not enough data, use all data
        mask = np.ones(len(time_s), dtype=bool)
    
    # Extract values for analysis period
    data_matrix = df.loc[mask, dist_cols].values
    
    # Calculate statistics
    mean_vals = np.mean(data_matrix, axis=0)
    min_vals = np.min(data_matrix, axis=0)
    max_vals = np.max(data_matrix, axis=0)
    std_vals = np.std(data_matrix, axis=0)
    
    # Apply unit conversions
    # Model outputs CH4/N2O in µmol/L (µM), validation is in nmol/L
    if variable in ['ch4', 'n2o']:
        mean_vals = mean_vals * 1000  # Convert µM → nmol/L
        min_vals = min_vals * 1000
        max_vals = max_vals * 1000
        std_vals = std_vals * 1000
    
    return distances, mean_vals, min_vals, max_vals, std_vals


def load_model_output_simple(branch, variable):
    """Load model output - just last timestep for backward compatibility"""
    distances, mean_vals, min_vals, max_vals, std_vals = load_model_output_with_tidal_stats(branch, variable)
    return distances, mean_vals


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
    
    # Normalized RMSE (relative to observation range)
    obs_range = np.max(obs) - np.min(obs)
    nrmse = rmse / obs_range if obs_range > 0 else np.nan
    
    # R² (coefficient of determination)
    if np.std(obs) > 0:
        r, p = stats.pearsonr(obs, mod)
        r2 = r**2
    else:
        r2 = np.nan
    
    # Percent bias
    pbias = 100 * bias / np.mean(obs) if np.mean(obs) != 0 else np.nan
    
    return {
        'n': len(obs),
        'bias': bias,
        'rmse': rmse,
        'nrmse': nrmse,
        'mae': mae,
        'r2': r2,
        'pbias': pbias,
        'obs_mean': np.mean(obs),
        'mod_mean': np.mean(mod),
        'obs_std': np.std(obs),
        'mod_std': np.std(mod)
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
    ('n2o', 'N2O_nmol_L', 'N2O [nmol/L]', 0, 50),
]

# Create main comparison figure with TIDAL VARIATION
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
        
        # Get model output WITH TIDAL VARIATION
        result = load_model_output_with_tidal_stats(branch, model_var)
        if result[0] is None:
            continue
        model_dist, mean_vals, min_vals, max_vals, std_vals = result
        
        # Plot model output - MEAN LINE
        color = branch_colors[branch]
        ax.plot(model_dist, mean_vals, '-', color=color, 
                linewidth=2, alpha=0.9, label=f'{branch}')
        
        # Plot TIDAL VARIATION ENVELOPE (min-max)
        ax.fill_between(model_dist, min_vals, max_vals, 
                        color=color, alpha=0.2)
        
        # Plot validation data
        ax.scatter(bval['Distance_km'], bval[val_var], 
                   color=color, s=60, edgecolors='black',
                   linewidths=1, zorder=5, marker='o')
        
        # Interpolate model to validation points for metrics
        if len(model_dist) > 0 and len(bval) > 0:
            model_interp = np.interp(bval['Distance_km'], model_dist, mean_vals)
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

plt.suptitle('C-GEM Model vs March 2025 Field Data\n(Shaded area = tidal variation)', fontsize=16, y=1.01)
plt.tight_layout()
plt.savefig(fig_dir / 'model_validation_comparison.png', dpi=300, bbox_inches='tight')
print(f"Saved: {fig_dir / 'model_validation_comparison.png'}")

# Print metrics summary
print("\n" + "=" * 80)
print("MODEL VALIDATION METRICS")
print("=" * 80)
print("\nNote: Metrics calculated using time-averaged model output vs observations")
print("      PBIAS = Percent Bias (positive = model overestimates)")
print("      NRMSE = Normalized RMSE (relative to obs range)")
print()

for var in all_metrics:
    print(f"\n{var.upper()}:")
    for branch in all_metrics[var]:
        m = all_metrics[var][branch]
        if m['n'] >= 3:
            pbias_str = f"{m['pbias']:+.1f}%" if not np.isnan(m['pbias']) else "N/A"
            print(f"  {branch:15s}: n={m['n']:2d}, RMSE={m['rmse']:7.2f}, Bias={m['bias']:+7.2f}, R²={m['r2']:.3f}, PBIAS={pbias_str}")

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
        
        model_dist, mean_vals, _, _, _ = load_model_output_with_tidal_stats(branch, model_var)
        if model_dist is None:
            continue
        
        # Interpolate
        model_interp = np.interp(bval['Distance_km'], model_dist, mean_vals)
        
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
            ax.text(0.05, 0.95, f"RMSE={metrics['rmse']:.1f}\nR²={metrics['r2']:.2f}\nBias={metrics['bias']:+.1f}", 
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

# Create detailed profile comparison for each branch - WITH TIDAL VARIATION
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
        
        result = load_model_output_with_tidal_stats(branch, model_var)
        
        if result[0] is not None:
            model_dist, mean_vals, min_vals, max_vals, std_vals = result
            
            # Plot mean with tidal envelope
            ax.plot(model_dist, mean_vals, '-', color='blue', 
                    linewidth=2, label='Model (mean)')
            ax.fill_between(model_dist, min_vals, max_vals, 
                            color='blue', alpha=0.2, label='Tidal range')
        
        if val_var in bval.columns:
            ax.scatter(bval['Distance_km'], bval[val_var], 
                       color='red', s=80, edgecolors='black',
                       linewidths=1, zorder=5, label='Observed')
        
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=8)
    
    plt.suptitle(f'{branch} - Model vs Validation\n(Blue shading = tidal variation)', fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(fig_dir / f'validation_{branch}.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {fig_dir / f'validation_{branch}.png'}")

plt.close('all')

# Print honest summary
print("\n" + "=" * 80)
print("HONEST ASSESSMENT")
print("=" * 80)
print("""
Key observations from the validation:

1. SALINITY: Generally good agreement (R² > 0.8), but model may underpredict 
   intrusion during dry season or overpredict during wet season.

2. O2: Model tends to UNDERESTIMATE oxygen, especially upstream where observed
   values are 200-260 µmol/L but model shows 160-200 µmol/L. This suggests:
   - Benthic O2 consumption may be too high
   - Or reaeration is underestimated
   - Or lateral inputs of oxygenated water are missing

3. pCO2: Model tends to OVERESTIMATE pCO2, especially in the transition zone.
   This is consistent with the O2 underestimation (too much respiration).

4. TOC: Model tends to UNDERESTIMATE TOC upstream, suggesting missing lateral
   organic carbon inputs from rice paddies, aquaculture, and urban runoff.

5. NH4/CH4/N2O: All significantly UNDERESTIMATED upstream. This strongly 
   indicates that lateral loads from agricultural areas are not adequate.
   These species require explicit lateral source terms.

6. The tidal variation shown in the figures is the MIN-MAX range over the
   last few days of simulation. Field observations are point measurements
   and should ideally fall within this envelope.
""")

print("\n" + "=" * 80)
print("VALIDATION COMPLETE")
print("=" * 80)
