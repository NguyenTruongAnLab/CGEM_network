#!/usr/bin/env python3
"""
Deep Analysis of Mekong Delta Field Data (March 2025)
======================================================

This script:
1. Deeply analyzes the validation data to understand estuarine patterns
2. Converts all units to match C-GEM model units
3. Creates a properly formatted CSV for model comparison
4. Generates diagnostic figures to understand the system

Author: C-GEM Team
Date: December 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE_DIR = Path('INPUT/Cases/Mekong_Delta_Full')
VAL_DIR = BASE_DIR / 'Validation'
OUTPUT_DIR = Path('OUTPUT/Mekong_Delta_Full/figures')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# UNIT CONVERSION FACTORS
# =============================================================================
# Field data units -> Model units (C-GEM)
CONVERSIONS = {
    # Dissolved gases
    'O2_mg_to_umol': 31.25,       # mg/L O2 -> µmol/L (MW O2 = 32 g/mol)
    
    # Nutrients (as N or P)
    'NO3N_mg_to_umol': 71.43,     # mgN/L -> µmol/L (MW N = 14 g/mol)
    'NH4N_mg_to_umol': 71.43,     # mgN/L -> µmol/L
    'PO4P_mg_to_umol': 32.26,     # mgP/L -> µmol/L (MW P = 31 g/mol)
    
    # Silicon
    'DSi_mg_to_umol': 35.61,      # mgSi/L -> µmol/L (MW Si = 28.09 g/mol)
    
    # Organic carbon
    'DOC_mg_to_umol': 83.33,      # mgC/L -> µmol/L (MW C = 12 g/mol)
    
    # Alkalinity
    'Alk_mmol_to_umol': 1000.0,   # mmol/L -> µmol/L
    
    # Chlorophyll-a to phytoplankton carbon
    # C:Chl ratio typically 40-60 for diatoms, 50-80 for others
    'Chla_ug_to_umolC': 0.0333,   # µg Chl-a/L -> µmol C/L (using C:Chl = 40)
    
    # Greenhouse gases
    # CH4: µgC/L -> nmol/L: divide by 12 (MW C), multiply by 1000
    'CH4_ugC_to_nmol': 83.33,     # µgC-CH4/L -> nmol/L
    
    # N2O: µgN/L -> nmol/L: divide by 28 (MW N2 in N2O), multiply by 1000
    'N2O_ugN_to_nmol': 35.71,     # µgN-N2O/L -> nmol/L
}

# =============================================================================
# Load and analyze data
# =============================================================================
print("=" * 80)
print("DEEP ANALYSIS OF MEKONG DELTA FIELD DATA (March 2025)")
print("=" * 80)

# Load raw data
df_raw = pd.read_csv(VAL_DIR / 'Mekong_4_branches_Mar_2025.csv')
print(f"\nLoaded {len(df_raw)} data points from 4 branches\n")

# Show raw data structure
print("RAW DATA COLUMNS:")
print("-" * 40)
for col in df_raw.columns:
    non_null = df_raw[col].notna().sum()
    if df_raw[col].dtype in ['float64', 'int64']:
        min_val = df_raw[col].min()
        max_val = df_raw[col].max()
        print(f"  {col:25s}: {non_null:3d} values, range [{min_val:.3f} - {max_val:.3f}]")
    else:
        print(f"  {col:25s}: {non_null:3d} values")

# =============================================================================
# Analyze each branch separately
# =============================================================================
print("\n" + "=" * 80)
print("BRANCH-BY-BRANCH ANALYSIS")
print("=" * 80)

branches = ['Hau_River', 'Ham_Luong', 'Co_Chien', 'My_Tho']
branch_stats = {}

for branch in branches:
    bdf = df_raw[df_raw['River'] == branch].sort_values('Distance')
    branch_stats[branch] = {
        'n_points': len(bdf),
        'dist_range': (bdf['Distance'].min(), bdf['Distance'].max()),
    }
    
    print(f"\n{'='*60}")
    print(f"BRANCH: {branch}")
    print(f"{'='*60}")
    print(f"  Points: {len(bdf)}")
    print(f"  Distance range: {bdf['Distance'].min():.1f} - {bdf['Distance'].max():.1f} km")
    
    # Salinity profile
    print(f"\n  SALINITY PROFILE (Surface Salinity, PSU):")
    for _, row in bdf.iterrows():
        sal = row['Surface Salinity']
        if pd.notna(sal):
            print(f"    {row['Distance']:6.1f} km: {sal:5.1f} PSU")
    
    # Find salinity intrusion length (where S < 1 PSU)
    fresh_zone = bdf[bdf['Surface Salinity'] < 1.0]
    if len(fresh_zone) > 0:
        intrusion_km = fresh_zone['Distance'].min()
        print(f"  → Salinity intrusion limit (< 1 PSU): ~{intrusion_km:.0f} km from mouth")
        branch_stats[branch]['intrusion_km'] = intrusion_km
    
    # Oxygen profile
    print(f"\n  OXYGEN PROFILE:")
    print(f"    {'Dist':>6s} {'O2(mg/L)':>10s} {'O2(%)':>8s} {'O2(µmol/L)':>12s}")
    for _, row in bdf.iterrows():
        o2_mg = row['O2 (mg/L)']
        o2_pct = row['O2 (%)']
        o2_umol = o2_mg * CONVERSIONS['O2_mg_to_umol'] if pd.notna(o2_mg) else np.nan
        if pd.notna(o2_mg):
            print(f"    {row['Distance']:6.1f} {o2_mg:10.2f} {o2_pct:8.1f} {o2_umol:12.1f}")
    
    # Find O2 minimum
    o2_min_idx = bdf['O2 (mg/L)'].idxmin()
    if pd.notna(o2_min_idx):
        o2_min = bdf.loc[o2_min_idx]
        print(f"  → O2 minimum: {o2_min['O2 (mg/L)']:.2f} mg/L ({o2_min['O2 (mg/L)']*31.25:.0f} µmol/L) at {o2_min['Distance']:.0f} km")
        branch_stats[branch]['o2_min_umol'] = o2_min['O2 (mg/L)'] * 31.25
        branch_stats[branch]['o2_min_km'] = o2_min['Distance']

    # pCO2 profile
    print(f"\n  pCO2 PROFILE:")
    print(f"    {'Dist':>6s} {'pCO2(ppm)':>12s} {'pH':>6s}")
    for _, row in bdf.iterrows():
        pco2 = row['pCO2 (ppm)']
        ph = row['pH']
        if pd.notna(pco2):
            print(f"    {row['Distance']:6.1f} {pco2:12.0f} {ph:6.2f}")
    
    # Find pCO2 maximum
    pco2_max_idx = bdf['pCO2 (ppm)'].idxmax()
    if pd.notna(pco2_max_idx):
        pco2_max = bdf.loc[pco2_max_idx]
        print(f"  → pCO2 maximum: {pco2_max['pCO2 (ppm)']:.0f} ppm at {pco2_max['Distance']:.0f} km (pH={pco2_max['pH']:.2f})")
        branch_stats[branch]['pco2_max'] = pco2_max['pCO2 (ppm)']
        branch_stats[branch]['pco2_max_km'] = pco2_max['Distance']

# =============================================================================
# Summary Statistics
# =============================================================================
print("\n" + "=" * 80)
print("OVERALL SUMMARY STATISTICS")
print("=" * 80)

# Ocean endmember (at mouth, distance ~ 0)
ocean_points = df_raw[df_raw['Distance'] < 1.0]
print("\n1. OCEAN ENDMEMBER (distance < 1 km from mouth):")
print(f"   Salinity: {ocean_points['Surface Salinity'].mean():.1f} ± {ocean_points['Surface Salinity'].std():.1f} PSU")
print(f"   O2:       {ocean_points['O2 (mg/L)'].mean():.2f} ± {ocean_points['O2 (mg/L)'].std():.2f} mg/L = {ocean_points['O2 (mg/L)'].mean()*31.25:.0f} µmol/L")
print(f"   pH:       {ocean_points['pH'].mean():.2f} ± {ocean_points['pH'].std():.2f}")
print(f"   pCO2:     {ocean_points['pCO2 (ppm)'].mean():.0f} ± {ocean_points['pCO2 (ppm)'].std():.0f} ppm")
print(f"   Alkalinity: {ocean_points['Alkalinity (mmol/L)'].mean():.2f} mmol/L = {ocean_points['Alkalinity (mmol/L)'].mean()*1000:.0f} µmol/L")

# River endmember (at upstream, distance > 60 km, salinity < 0.5)
river_points = df_raw[(df_raw['Distance'] > 60) & (df_raw['Surface Salinity'] < 0.5)]
print("\n2. RIVER ENDMEMBER (distance > 60 km, salinity < 0.5 PSU):")
print(f"   N points: {len(river_points)}")
print(f"   Salinity: {river_points['Surface Salinity'].mean():.2f} PSU (essentially fresh)")
print(f"   O2:       {river_points['O2 (mg/L)'].mean():.2f} ± {river_points['O2 (mg/L)'].std():.2f} mg/L = {river_points['O2 (mg/L)'].mean()*31.25:.0f} µmol/L")
print(f"   pH:       {river_points['pH'].mean():.2f} ± {river_points['pH'].std():.2f}")
print(f"   pCO2:     {river_points['pCO2 (ppm)'].mean():.0f} ± {river_points['pCO2 (ppm)'].std():.0f} ppm")
print(f"   Alkalinity: {river_points['Alkalinity (mmol/L)'].mean():.2f} mmol/L = {river_points['Alkalinity (mmol/L)'].mean()*1000:.0f} µmol/L")
print(f"   NO3:      {river_points['NO3-N (mgN/L)'].mean():.2f} mgN/L = {river_points['NO3-N (mgN/L)'].mean()*71.43:.1f} µmol/L")
print(f"   NH4:      {river_points['NH4 (mgN/L)'].mean():.3f} mgN/L = {river_points['NH4 (mgN/L)'].mean()*71.43:.1f} µmol/L")
print(f"   DOC:      {river_points['DOC (mgC/L)'].mean():.2f} mgC/L = {river_points['DOC (mgC/L)'].mean()*83.33:.0f} µmol/L")
print(f"   TSS:      {river_points['TSS (mg/L)'].mean():.1f} mg/L")
print(f"   CH4:      {river_points['ugC-CH4/L'].mean():.2f} µgC/L = {river_points['ugC-CH4/L'].mean()*83.33:.0f} nmol/L")
print(f"   N2O:      {river_points['ugN-N2O/L'].mean():.2f} µgN/L = {river_points['ugN-N2O/L'].mean()*35.71:.0f} nmol/L")

# Mid-estuary (transition zone)
mid_points = df_raw[(df_raw['Distance'] > 20) & (df_raw['Distance'] < 50)]
print("\n3. MID-ESTUARY (20-50 km from mouth):")
print(f"   N points: {len(mid_points)}")
print(f"   Salinity: {mid_points['Surface Salinity'].mean():.1f} ± {mid_points['Surface Salinity'].std():.1f} PSU")
print(f"   O2:       {mid_points['O2 (mg/L)'].mean():.2f} ± {mid_points['O2 (mg/L)'].std():.2f} mg/L")
print(f"   TSS:      {mid_points['TSS (mg/L)'].mean():.1f} ± {mid_points['TSS (mg/L)'].std():.1f} mg/L (ETM zone?)")

# =============================================================================
# Create converted CSV file with model units
# =============================================================================
print("\n" + "=" * 80)
print("CREATING CONVERTED CSV FILE")
print("=" * 80)

# Create new dataframe with model units
df_model = pd.DataFrame()

# Keep identifiers
df_model['Station'] = df_raw['Name']
df_model['Branch'] = df_raw['River']
df_model['Distance_km'] = df_raw['Distance']
df_model['Year'] = df_raw['Year']

# Temperature (already in °C)
df_model['Temperature_C'] = df_raw['T°']

# Salinity (already in PSU)
df_model['Salinity_PSU'] = df_raw['Surface Salinity']

# pH (dimensionless)
df_model['pH'] = df_raw['pH']

# Oxygen (mg/L -> µmol/L)
df_model['O2_umol_L'] = df_raw['O2 (mg/L)'] * CONVERSIONS['O2_mg_to_umol']
df_model['O2_percent'] = df_raw['O2 (%)']

# TSS/SPM (already in mg/L - model uses mg/L)
df_model['SPM_mg_L'] = df_raw['TSS (mg/L)']

# Nutrients
df_model['NO3_umol_L'] = df_raw['NO3-N (mgN/L)'] * CONVERSIONS['NO3N_mg_to_umol']
df_model['NH4_umol_L'] = df_raw['NH4 (mgN/L)'] * CONVERSIONS['NH4N_mg_to_umol']
df_model['PO4_umol_L'] = df_raw['PO4 (mgP/L)'] * CONVERSIONS['PO4P_mg_to_umol']
df_model['DSi_umol_L'] = df_raw['Dsi (mgSi/L)'] * CONVERSIONS['DSi_mg_to_umol']

# Organic carbon (DOC to TOC proxy)
df_model['TOC_umol_L'] = df_raw['DOC (mgC/L)'] * CONVERSIONS['DOC_mg_to_umol']

# Chlorophyll-a -> Phytoplankton carbon (approximation)
df_model['Chla_ug_L'] = df_raw['Chl-a (ug/L)']
df_model['PHY_umolC_L'] = df_raw['Chl-a (ug/L)'] * CONVERSIONS['Chla_ug_to_umolC']

# Alkalinity (mmol/L -> µmol/L)
df_model['Alkalinity_umol_L'] = df_raw['Alkalinity (mmol/L)'] * CONVERSIONS['Alk_mmol_to_umol']

# pCO2 (already in ppm)
df_model['pCO2_ppm'] = df_raw['pCO2 (ppm)']

# Greenhouse gases
df_model['CH4_nmol_L'] = df_raw['ugC-CH4/L'] * CONVERSIONS['CH4_ugC_to_nmol']
df_model['N2O_nmol_L'] = df_raw['ugN-N2O/L'] * CONVERSIONS['N2O_ugN_to_nmol']

# Save converted file
output_csv = VAL_DIR / 'Mekong_Mar2025_model_units.csv'
df_model.to_csv(output_csv, index=False, float_format='%.2f')
print(f"\nSaved converted data to: {output_csv}")

# Print converted data summary
print("\nCONVERTED DATA SUMMARY (Model Units):")
print("-" * 60)
for col in df_model.columns[4:]:  # Skip identifiers
    non_null = df_model[col].notna().sum()
    if non_null > 0:
        print(f"  {col:25s}: min={df_model[col].min():8.2f}, max={df_model[col].max():8.2f}, mean={df_model[col].mean():8.2f}")

# =============================================================================
# Create summary statistics for model boundary conditions
# =============================================================================
print("\n" + "=" * 80)
print("RECOMMENDED MODEL BOUNDARY CONDITIONS")
print("=" * 80)

print("\n1. OCEAN BOUNDARY (downstream, all tidal outlets):")
print("   [Based on mouth samples, distance < 2 km]")
ocean = df_model[df_model['Distance_km'] < 2]
print(f"   Salinity:    {ocean['Salinity_PSU'].mean():.1f} PSU")
print(f"   O2:          {ocean['O2_umol_L'].mean():.0f} µmol/L")
print(f"   pH:          {ocean['pH'].mean():.2f}")
print(f"   pCO2:        {ocean['pCO2_ppm'].mean():.0f} ppm")
print(f"   Alkalinity:  {ocean['Alkalinity_umol_L'].mean():.0f} µmol/L")
print(f"   NO3:         {ocean['NO3_umol_L'].mean():.1f} µmol/L")
print(f"   NH4:         {ocean['NH4_umol_L'].mean():.1f} µmol/L")
print(f"   PO4:         {ocean['PO4_umol_L'].mean():.1f} µmol/L")
print(f"   SPM:         {ocean['SPM_mg_L'].mean():.0f} mg/L")
print(f"   TOC:         {ocean['TOC_umol_L'].mean():.0f} µmol/L")
print(f"   CH4:         {ocean['CH4_nmol_L'].mean():.0f} nmol/L")
print(f"   N2O:         {ocean['N2O_nmol_L'].mean():.0f} nmol/L")

print("\n2. RIVER BOUNDARY (upstream, discharge nodes):")
print("   [Based on freshwater samples, distance > 65 km, salinity < 0.5 PSU]")
river = df_model[(df_model['Distance_km'] > 65) & (df_model['Salinity_PSU'] < 0.5)]
print(f"   Salinity:    {river['Salinity_PSU'].mean():.2f} PSU")
print(f"   O2:          {river['O2_umol_L'].mean():.0f} µmol/L")
print(f"   pH:          {river['pH'].mean():.2f}")
print(f"   pCO2:        {river['pCO2_ppm'].mean():.0f} ppm")
print(f"   Alkalinity:  {river['Alkalinity_umol_L'].mean():.0f} µmol/L")
print(f"   NO3:         {river['NO3_umol_L'].mean():.1f} µmol/L")
print(f"   NH4:         {river['NH4_umol_L'].mean():.1f} µmol/L")
print(f"   PO4:         {river['PO4_umol_L'].mean():.1f} µmol/L")
print(f"   SPM:         {river['SPM_mg_L'].mean():.0f} mg/L")
print(f"   TOC:         {river['TOC_umol_L'].mean():.0f} µmol/L")
print(f"   CH4:         {river['CH4_nmol_L'].mean():.0f} nmol/L")
print(f"   N2O:         {river['N2O_nmol_L'].mean():.0f} nmol/L")

# =============================================================================
# Create diagnostic figures
# =============================================================================
print("\n" + "=" * 80)
print("CREATING DIAGNOSTIC FIGURES")
print("=" * 80)

# Figure 1: Salinity profiles for all branches
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

colors = {'Hau_River': '#ff7f0e', 'Ham_Luong': '#e377c2', 'My_Tho': '#8c564b', 'Co_Chien': '#d62728'}

for idx, branch in enumerate(branches):
    ax = axes[idx]
    bdf = df_model[df_model['Branch'] == branch].sort_values('Distance_km')
    
    # Primary: Salinity
    ax.plot(bdf['Distance_km'], bdf['Salinity_PSU'], 'o-', color=colors[branch], 
            linewidth=2, markersize=8, label='Salinity')
    ax.axhline(y=4, color='red', linestyle='--', alpha=0.5, label='4 PSU threshold')
    ax.axhline(y=1, color='orange', linestyle='--', alpha=0.5, label='1 PSU threshold')
    
    ax.set_xlabel('Distance from mouth (km)')
    ax.set_ylabel('Salinity (PSU)', color=colors[branch])
    ax.set_title(f'{branch}')
    ax.set_xlim(0, bdf['Distance_km'].max() + 5)
    ax.set_ylim(0, 35)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

plt.suptitle('Salinity Profiles - Mekong Delta (March 2025)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'validation_salinity_profiles.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'validation_salinity_profiles.png'}")

# Figure 2: O2 and pCO2 profiles
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for idx, branch in enumerate(branches):
    ax = axes[idx]
    bdf = df_model[df_model['Branch'] == branch].sort_values('Distance_km')
    
    # O2 (left axis)
    ax.plot(bdf['Distance_km'], bdf['O2_umol_L'], 'o-', color='blue', 
            linewidth=2, markersize=8, label='O2')
    ax.set_xlabel('Distance from mouth (km)')
    ax.set_ylabel('O2 (µmol/L)', color='blue')
    ax.tick_params(axis='y', labelcolor='blue')
    ax.set_ylim(100, 300)
    
    # pCO2 (right axis)
    ax2 = ax.twinx()
    ax2.plot(bdf['Distance_km'], bdf['pCO2_ppm'], 's-', color='red', 
             linewidth=2, markersize=6, label='pCO2')
    ax2.set_ylabel('pCO2 (ppm)', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_ylim(0, 5500)
    ax2.axhline(y=420, color='green', linestyle='--', alpha=0.5)
    
    ax.set_title(f'{branch}')
    ax.set_xlim(0, bdf['Distance_km'].max() + 5)
    ax.grid(True, alpha=0.3)
    
    # Combined legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2 + ['Atm CO2'], loc='upper right')

plt.suptitle('O2 and pCO2 Profiles - Mekong Delta (March 2025)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'validation_o2_pco2_profiles.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'validation_o2_pco2_profiles.png'}")

# Figure 3: SPM and Chl-a profiles
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for idx, branch in enumerate(branches):
    ax = axes[idx]
    bdf = df_model[df_model['Branch'] == branch].sort_values('Distance_km')
    
    # SPM (left axis)
    ax.fill_between(bdf['Distance_km'], 0, bdf['SPM_mg_L'], color='brown', alpha=0.3)
    ax.plot(bdf['Distance_km'], bdf['SPM_mg_L'], 'o-', color='brown', 
            linewidth=2, markersize=8, label='SPM')
    ax.set_xlabel('Distance from mouth (km)')
    ax.set_ylabel('SPM (mg/L)', color='brown')
    ax.tick_params(axis='y', labelcolor='brown')
    ax.set_ylim(0, 45)
    
    # Chl-a (right axis)
    ax2 = ax.twinx()
    ax2.plot(bdf['Distance_km'], bdf['Chla_ug_L'], 's-', color='green', 
             linewidth=2, markersize=6, label='Chl-a')
    ax2.set_ylabel('Chl-a (µg/L)', color='green')
    ax2.tick_params(axis='y', labelcolor='green')
    ax2.set_ylim(0, 8)
    
    ax.set_title(f'{branch}')
    ax.set_xlim(0, bdf['Distance_km'].max() + 5)
    ax.grid(True, alpha=0.3)
    
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.suptitle('SPM and Chlorophyll-a Profiles - Mekong Delta (March 2025)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'validation_spm_chla_profiles.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'validation_spm_chla_profiles.png'}")

# Figure 4: Nutrient profiles
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for idx, branch in enumerate(branches):
    ax = axes[idx]
    bdf = df_model[df_model['Branch'] == branch].sort_values('Distance_km')
    
    ax.plot(bdf['Distance_km'], bdf['NO3_umol_L'], 'o-', color='blue', 
            linewidth=2, markersize=6, label='NO3')
    ax.plot(bdf['Distance_km'], bdf['NH4_umol_L'], 's-', color='orange', 
            linewidth=2, markersize=6, label='NH4')
    ax.plot(bdf['Distance_km'], bdf['PO4_umol_L'], '^-', color='green', 
            linewidth=2, markersize=6, label='PO4')
    ax.plot(bdf['Distance_km'], bdf['DSi_umol_L'], 'd-', color='purple', 
            linewidth=2, markersize=6, label='DSi')
    
    ax.set_xlabel('Distance from mouth (km)')
    ax.set_ylabel('Concentration (µmol/L)')
    ax.set_title(f'{branch}')
    ax.set_xlim(0, bdf['Distance_km'].max() + 5)
    ax.set_ylim(0, 80)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

plt.suptitle('Nutrient Profiles - Mekong Delta (March 2025)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'validation_nutrient_profiles.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'validation_nutrient_profiles.png'}")

# Figure 5: GHG profiles
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for idx, branch in enumerate(branches):
    ax = axes[idx]
    bdf = df_model[df_model['Branch'] == branch].sort_values('Distance_km')
    
    # CH4 (left axis)
    ax.plot(bdf['Distance_km'], bdf['CH4_nmol_L'], 'o-', color='orange', 
            linewidth=2, markersize=8, label='CH4')
    ax.set_xlabel('Distance from mouth (km)')
    ax.set_ylabel('CH4 (nmol/L)', color='orange')
    ax.tick_params(axis='y', labelcolor='orange')
    ax.set_ylim(0, 300)
    
    # N2O (right axis)
    ax2 = ax.twinx()
    ax2.plot(bdf['Distance_km'], bdf['N2O_nmol_L'], 's-', color='purple', 
             linewidth=2, markersize=6, label='N2O')
    ax2.set_ylabel('N2O (nmol/L)', color='purple')
    ax2.tick_params(axis='y', labelcolor='purple')
    ax2.set_ylim(0, 200)
    
    ax.set_title(f'{branch}')
    ax.set_xlim(0, bdf['Distance_km'].max() + 5)
    ax.grid(True, alpha=0.3)
    
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.suptitle('Greenhouse Gas Profiles - Mekong Delta (March 2025)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'validation_ghg_profiles.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'validation_ghg_profiles.png'}")

# Figure 6: pH and Alkalinity
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for idx, branch in enumerate(branches):
    ax = axes[idx]
    bdf = df_model[df_model['Branch'] == branch].sort_values('Distance_km')
    
    # pH (left axis)
    ax.plot(bdf['Distance_km'], bdf['pH'], 'o-', color='blue', 
            linewidth=2, markersize=8, label='pH')
    ax.set_xlabel('Distance from mouth (km)')
    ax.set_ylabel('pH', color='blue')
    ax.tick_params(axis='y', labelcolor='blue')
    ax.set_ylim(7.2, 8.3)
    
    # Alkalinity (right axis)
    ax2 = ax.twinx()
    ax2.plot(bdf['Distance_km'], bdf['Alkalinity_umol_L'], 's-', color='red', 
             linewidth=2, markersize=6, label='Alkalinity')
    ax2.set_ylabel('Alkalinity (µmol/L)', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_ylim(1200, 2400)
    
    ax.set_title(f'{branch}')
    ax.set_xlim(0, bdf['Distance_km'].max() + 5)
    ax.grid(True, alpha=0.3)
    
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.suptitle('pH and Alkalinity Profiles - Mekong Delta (March 2025)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'validation_ph_alk_profiles.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'validation_ph_alk_profiles.png'}")

plt.close('all')

# =============================================================================
# KEY INSIGHTS FOR MODEL
# =============================================================================
print("\n" + "=" * 80)
print("KEY INSIGHTS FOR C-GEM MODEL")
print("=" * 80)

print("""
1. SALINITY INTRUSION:
   - Sharp gradient in mid-estuary (10-40 km)
   - Intrusion lengths: Hau ~42 km, Ham Luong ~50 km, Co Chien ~55 km, My Tho ~36 km
   - Ocean endmember: ~30 PSU

2. OXYGEN DYNAMICS:
   - Ocean: ~260 µmol/L (saturated, ~105%)
   - River: ~170-180 µmol/L (undersaturated, ~70-75%)
   - Clear DECREASE upstream → indicates net respiration exceeds reaeration
   - This is OPPOSITE to what the model currently shows!

3. pCO2 DYNAMICS:
   - Ocean: ~500 ppm (near atmospheric)
   - River: ~4000-4700 ppm (HIGHLY supersaturated)
   - Clear INCREASE upstream → consistent with respiration
   - pCO2 and O2 are inversely correlated (as expected)

4. NUTRIENTS:
   - Very LOW throughout (oligotrophic conditions)
   - NO3: 6-70 µmol/L (highest at ocean, diluted upstream!)
   - NH4: 0.7-5 µmol/L (very low)
   - This suggests: (1) low nutrient river, (2) biological uptake

5. SPM/TSS:
   - 6-38 mg/L (moderate turbidity)
   - ETM zone visible in some branches (20-40 km)

6. CARBONATE CHEMISTRY:
   - pH: 8.1 (ocean) → 7.4 (river) - consistent with pCO2 increase
   - Alkalinity: ~2100 µmol/L (ocean) → ~1300 µmol/L (river)
   - This explains the pCO2 pattern: lower buffering capacity upstream

7. GHGs:
   - CH4: 40-240 nmol/L (increases upstream - benthic sources)
   - N2O: 7-190 nmol/L (variable, some high values near Can Tho)
""")

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
