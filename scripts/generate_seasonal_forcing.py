#!/usr/bin/env python3
"""
Generate Seasonal Forcing Data for Mekong Delta Full Case
==========================================================

This script generates SEASONALLY-VARYING forcing data files that capture
the critical dry-wet cycle of the Mekong Delta:

MEKONG SEASONAL HYDROLOGY (MRC 2018, Nguyen et al. 2008):
=========================================================
DRY SEASON (December - May):
- Low discharge: ~2,000-4,000 m³/s at Vietnamese border
- Maximum salt intrusion: 40-60 km inland
- Low SPM: 20-50 mg/L (clear water)
- Higher nutrient concentrations (less dilution)

WET SEASON (June - November):
- High discharge: 15,000-40,000 m³/s (flood peak in September)
- Minimal salt intrusion: <10 km
- High SPM: 150-400 mg/L (monsoon runoff)
- Lower nutrient concentrations (dilution)

SPECIES VARIATIONS:
- SPM: 5x higher in wet season (monsoon sediment load)
- NO3/NH4: 1.5x higher in dry season (concentration effect)
- DO: Slightly lower in wet season (higher respiration from organic load)
- Salinity: Ocean BC stays constant; river BC always ~0.1 PSU

Output Directory: INPUT/Cases/Mekong_Delta_Full/forcing_data/

Author: Nguyen Truong An
Date: December 2025
"""

import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta

# ===========================================================================
# CONFIGURATION
# ===========================================================================

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent
CASE_DIR = PROJECT_ROOT / "INPUT" / "Cases" / "Mekong_Delta_Full"
FORCING_DIR = CASE_DIR / "forcing_data"

# Simulation parameters (must match case_config.txt)
START_DATE = "2024-01-01"
DURATION_DAYS = 465  # ~15 months to capture full seasonal cycle
DT_SPECIES = 86400.0  # Daily species [s]

# ===========================================================================
# SEASONAL PARAMETERS
# ===========================================================================

# Month-specific discharge multipliers (relative to annual mean)
# Based on MRC Kratie gauge data
MONTHLY_Q_MULTIPLIER = {
    1: 0.35,   # January - Dry
    2: 0.28,   # February - Driest
    3: 0.25,   # March - Driest
    4: 0.30,   # April - Late dry
    5: 0.50,   # May - Transition
    6: 0.90,   # June - Early wet
    7: 1.40,   # July - Wet
    8: 2.00,   # August - Flood rising
    9: 2.50,   # September - Flood peak
    10: 2.20,  # October - Flood falling
    11: 1.20,  # November - Late wet
    12: 0.55,  # December - Transition to dry
}

# Annual mean discharge at Vietnamese border [m³/s]
# Standardized to match setup_mekong_simplified.py (based on MRC Kratie average)
Q_ANNUAL_MEAN = 13000.0

# ===========================================================================
# RIVER SPECIES CONCENTRATIONS (Seasonally Varying)
# Based on Li et al. (2014), MRC monitoring, Hung & Hieu (2019)
# ===========================================================================

# DRY SEASON values (concentrated)
SPECIES_RIVER_DRY = {
    "salinity": 0.1,       # PSU (always fresh)
    "phy1": 25.0,          # µgC/L (more phytoplankton in clear water)
    "phy2": 15.0,          # µgC/L
    "dsi": 120.0,          # µmol/L (higher in dry - basalt weathering)
    "no3": 50.0,           # µmolN/L (concentrated)
    "nh4": 18.0,           # µmolN/L (concentrated)
    "po4": 1.2,            # µmolP/L (moderate - realistic for Mekong)
    "o2": 250.0,           # µmol/L (well-oxygenated)
    "toc": 250.0,          # µmolC/L (moderate - less terrestrial input)
    "spm": 30.0,           # mg/L (clear water)
    "dic": 1950.0,         # µmol/L
    "at": 2050.0,          # µeq/L
    "pco2": 0.0,           # computed
    "co2": 0.0,            # computed
    "ph": 7.6,             # pH units
    "hs": 0.3,             # µmol/L
}

# WET SEASON values (diluted, high sediment)
SPECIES_RIVER_WET = {
    "salinity": 0.1,       # PSU (always fresh)
    "phy1": 10.0,          # µgC/L (light-limited by turbidity)
    "phy2": 5.0,           # µgC/L
    "dsi": 80.0,           # µmol/L (diluted)
    "no3": 30.0,           # µmolN/L (diluted)
    "nh4": 8.0,            # µmolN/L (diluted)
    "po4": 0.6,            # µmolP/L (diluted)
    "o2": 220.0,           # µmol/L (slightly lower due to organic load)
    "toc": 350.0,          # µmolC/L (higher - terrestrial input)
    "spm": 200.0,          # mg/L (monsoon sediment load)
    "dic": 1700.0,         # µmol/L (diluted)
    "at": 1800.0,          # µeq/L (diluted)
    "pco2": 0.0,           # computed
    "co2": 0.0,            # computed
    "ph": 7.3,             # pH units (more acidic with organic load)
    "hs": 0.5,             # µmol/L
}

# OCEAN SPECIES (constant - South China Sea)
SPECIES_OCEAN = {
    "salinity": 33.0,
    "phy1": 5.0,
    "phy2": 3.0,
    "dsi": 4.0,
    "no3": 1.5,
    "nh4": 0.3,
    "po4": 0.15,
    "o2": 210.0,
    "toc": 80.0,
    "spm": 8.0,
    "dic": 2100.0,
    "at": 2350.0,
    "pco2": 0.0,
    "co2": 0.0,
    "ph": 8.1,
    "hs": 0.1,
}


def get_seasonal_factor(day_of_year: int) -> float:
    """
    Get seasonal interpolation factor (0=dry, 1=wet).
    
    Uses sinusoidal approximation with peak wet in September (day 270).
    """
    # Peak wet at day 270 (late September)
    phase = 2 * np.pi * (day_of_year - 270) / 365.0
    factor = 0.5 * (1.0 - np.cos(phase))
    return factor


def interpolate_seasonal(dry_val: float, wet_val: float, factor: float) -> float:
    """Linearly interpolate between dry and wet season values."""
    return dry_val + factor * (wet_val - dry_val)


def generate_river_species_timeseries(duration_days: int) -> pd.DataFrame:
    """
    Generate seasonally-varying river species concentrations.
    """
    n_steps = int(duration_days * 86400 / DT_SPECIES) + 1
    time_s = np.arange(n_steps) * DT_SPECIES
    
    data = {"time_s": time_s}
    
    np.random.seed(42)
    
    for name in SPECIES_RIVER_DRY.keys():
        dry_val = SPECIES_RIVER_DRY[name]
        wet_val = SPECIES_RIVER_WET[name]
        
        values = []
        for t in time_s:
            # Convert time to day of year
            day = int((t / 86400) % 365)
            factor = get_seasonal_factor(day)
            
            # Interpolate between dry and wet
            base_val = interpolate_seasonal(dry_val, wet_val, factor)
            
            # Add small noise (±5%)
            noise = 1.0 + 0.05 * np.random.randn()
            values.append(base_val * np.clip(noise, 0.9, 1.1))
        
        data[name] = values
    
    return pd.DataFrame(data)


def generate_ocean_species_timeseries(duration_days: int) -> pd.DataFrame:
    """
    Generate ocean species concentrations (nearly constant).
    """
    n_steps = int(duration_days * 86400 / DT_SPECIES) + 1
    time_s = np.arange(n_steps) * DT_SPECIES
    
    data = {"time_s": time_s}
    
    np.random.seed(43)
    
    for name, val in SPECIES_OCEAN.items():
        # Small noise (±2%) - ocean is more stable
        noise = 1.0 + 0.02 * np.random.randn(n_steps)
        data[name] = val * np.clip(noise, 0.95, 1.05)
    
    return pd.DataFrame(data)


def generate_discharge_timeseries(duration_days: int, river: str) -> pd.DataFrame:
    """
    Generate seasonally-varying discharge time series.
    """
    dt = 3600.0  # Hourly
    n_steps = int(duration_days * 86400 / dt) + 1
    time_s = np.arange(n_steps) * dt
    
    # Tien/Hau split
    if river == "Tien":
        fraction = 0.60
    else:
        fraction = 0.40
    
    np.random.seed(44 if river == "Tien" else 45)
    
    Q_values = []
    for t in time_s:
        # Get month (1-12)
        day_of_year = int((t / 86400) % 365)
        month = (day_of_year // 30) % 12 + 1
        
        # Monthly multiplier
        mult = MONTHLY_Q_MULTIPLIER[month]
        
        # Base discharge
        Q_base = Q_ANNUAL_MEAN * mult * fraction
        
        # Daily variation (±10%)
        daily_var = 0.1 * Q_base * np.sin(2 * np.pi * t / 86400)
        
        # Small noise (±3%)
        noise = Q_base * 0.03 * np.random.randn()
        
        Q = Q_base + daily_var + noise
        Q_values.append(max(Q, Q_base * 0.5))  # Ensure positive
    
    return pd.DataFrame({"time_s": time_s, "Q_m3s": Q_values})


def main():
    """Generate all seasonal forcing data files."""
    
    print("=" * 70)
    print("Generating SEASONAL Forcing Data: Mekong Delta Full")
    print("=" * 70)
    print(f"Duration: {DURATION_DAYS} days (~{DURATION_DAYS/30:.1f} months)")
    print(f"Output directory: {FORCING_DIR}")
    print()
    
    FORCING_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. DISCHARGE FILES
    print("1. Generating seasonal discharge time series...")
    
    Q_tien = generate_discharge_timeseries(DURATION_DAYS, "Tien")
    Q_tien.to_csv(FORCING_DIR / "Tien_Inlet.csv", index=False)
    print(f"   - Tien_Inlet.csv: Q_range = {Q_tien['Q_m3s'].min():.0f} - {Q_tien['Q_m3s'].max():.0f} m³/s")
    
    Q_hau = generate_discharge_timeseries(DURATION_DAYS, "Hau")
    Q_hau.to_csv(FORCING_DIR / "Hau_Inlet.csv", index=False)
    print(f"   - Hau_Inlet.csv: Q_range = {Q_hau['Q_m3s'].min():.0f} - {Q_hau['Q_m3s'].max():.0f} m³/s")
    
    # 2. SPECIES BOUNDARY FILES
    print("\n2. Generating seasonal species boundary conditions...")
    
    # River species (seasonal)
    species_river = generate_river_species_timeseries(DURATION_DAYS)
    species_river.to_csv(FORCING_DIR / "species_river_ts.csv", index=False)
    print(f"   - species_river_ts.csv: {len(species_river)} time points")
    print(f"     SPM range: {species_river['spm'].min():.0f} - {species_river['spm'].max():.0f} mg/L")
    print(f"     NO3 range: {species_river['no3'].min():.0f} - {species_river['no3'].max():.0f} µmol/L")
    
    # Ocean species (constant)
    species_ocean = generate_ocean_species_timeseries(DURATION_DAYS)
    species_ocean.to_csv(FORCING_DIR / "species_ocean_ts.csv", index=False)
    print(f"   - species_ocean_ts.csv: {len(species_ocean)} time points")
    print(f"     Salinity = {SPECIES_OCEAN['salinity']:.0f} PSU (constant)")
    
    # Also update the root species files
    species_river.to_csv(CASE_DIR / "species_river.csv", index=False)
    species_ocean.to_csv(CASE_DIR / "species_ocean.csv", index=False)
    
    print("\n" + "=" * 70)
    print("Seasonal Forcing Generation Complete!")
    print("=" * 70)
    print()
    print("Seasonal summary:")
    print(f"  DRY season (Feb-Apr): Q ~{Q_ANNUAL_MEAN * 0.28 * 0.6:.0f} m³/s Tien, SPM ~30 mg/L")
    print(f"  WET season (Aug-Oct): Q ~{Q_ANNUAL_MEAN * 2.0 * 0.6:.0f} m³/s Tien, SPM ~200 mg/L")
    print()
    print("Next: Run the model to see salt intrusion respond to seasonal discharge")
    
    return 0


if __name__ == "__main__":
    exit(main())
