#!/usr/bin/env python3
"""
Generate Synthetic Forcing Data for Mekong Delta Full Case
==========================================================

This script generates all time-varying forcing data files:
- Discharge time series for Tien and Hau inlets
- Tidal elevation time series for 5 river mouths
- Species boundary condition files with time_s column

The data is synthetic but based on realistic Mekong Delta parameters:
- Dry season (March) discharge: ~4000-5000 m³/s total
- M2 tidal amplitude: 1.5-2.5 m
- Phase lag between mouths: ~8-10° per 25 km

Output Directory: INPUT/Cases/Mekong_Delta_Full/forcing_data/

Author: Nguyen Truong An
Date: November 2024
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
OUTPUT_DIR = PROJECT_ROOT / "OUTPUT" / "Mekong_Delta_Full"

# Simulation parameters (must match case_config.txt)
START_DATE = "2024-03-01"
DURATION_DAYS = 10
WARMUP_DAYS = 3
TOTAL_DAYS = DURATION_DAYS + WARMUP_DAYS

# Time resolution
DT_DISCHARGE = 3600.0     # Hourly discharge [s]
DT_TIDE = 1800.0          # Half-hourly tides [s]
DT_SPECIES = 86400.0      # Daily species [s]

# ===========================================================================
# DISCHARGE SPLIT (Academically Correct)
# ===========================================================================
# At Phnom Penh: 80% Tien, 20% Hau
# After Vam Nao redistribution (~30% transfer): 
#   - Tien: 60% effective
#   - Hau: 40% effective
# These fractions represent discharge at Vietnamese border

TIEN_FRACTION = 0.60      # Tien receives 60% of total
HAU_FRACTION = 0.40       # Hau receives 40% of total

# ===========================================================================
# HYDROLOGICAL PARAMETERS (Dry Season - March)
# ===========================================================================

# Total Mekong discharge at Vietnamese border [m³/s]
# Dry season (March): ~4000-5500 m³/s
Q_TOTAL_MEAN = 4500.0     # Mean discharge [m³/s]
Q_TOTAL_AMPLITUDE = 500.0  # Daily variation amplitude [m³/s]

# ===========================================================================
# TIDAL PARAMETERS
# ===========================================================================

# M2 tide
M2_PERIOD_S = 12.42 * 3600  # M2 period [s]
M2_OMEGA = 2 * np.pi / M2_PERIOD_S

# Spring-neap modulation (Msf constituent)
MSF_PERIOD_S = 14.77 * 86400  # 14.77 days
MSF_OMEGA = 2 * np.pi / MSF_PERIOD_S
MSF_AMPLITUDE = 0.2  # ±20% modulation

# Tidal parameters for each mouth (from Nguyen et al. 2008)
# Amplitude [m], Phase [degrees] relative to My Tho
TIDAL_MOUTHS = {
    "MyTho": {"amplitude": 1.8, "phase_deg": 0.0},
    "HamLuong": {"amplitude": 1.9, "phase_deg": 8.0},
    "CoChien": {"amplitude": 2.0, "phase_deg": 15.0},
    "DinhAn": {"amplitude": 2.2, "phase_deg": 25.0},
    "TranDe": {"amplitude": 2.0, "phase_deg": 35.0},
    # Merged Hau outlet (average of Dinh An + Tran De)
    "Hau": {"amplitude": 2.1, "phase_deg": 30.0},
}

# ===========================================================================
# SPECIES CONCENTRATIONS
# ===========================================================================

# River boundary (freshwater, nutrient-rich)
SPECIES_RIVER = {
    "salinity": 0.1,       # PSU
    "phy1": 15.0,          # µgC/L
    "phy2": 8.0,           # µgC/L
    "dsi": 80.0,           # µmol/L
    "no3": 35.0,           # µmolN/L
    "nh4": 12.0,           # µmolN/L
    "po4": 1.5,            # µmolP/L
    "o2": 250.0,           # µmol/L (well-oxygenated)
    "toc": 300.0,          # µmolC/L
    "spm": 80.0,           # mg/L
    "dic": 1800.0,         # µmol/L
    "at": 1900.0,          # µeq/L
    "pco2": 0.0,           # computed
    "co2": 0.0,            # computed
    "ph": 7.5,             # pH units
    "hs": 0.5,             # µmol/L
}

# Ocean boundary (seawater, oligotrophic)
SPECIES_OCEAN = {
    "salinity": 32.0,
    "phy1": 5.0,
    "phy2": 3.0,
    "dsi": 5.0,
    "no3": 2.0,
    "nh4": 0.5,
    "po4": 0.2,
    "o2": 210.0,
    "toc": 100.0,
    "spm": 10.0,
    "dic": 2100.0,
    "at": 2350.0,
    "pco2": 0.0,
    "co2": 0.0,
    "ph": 8.1,
    "hs": 0.1,
}


# ===========================================================================
# DATA GENERATORS
# ===========================================================================

def generate_discharge(duration_days: float, river: str = "Tien") -> pd.DataFrame:
    """
    Generate discharge time series.
    
    Args:
        duration_days: Simulation duration [days]
        river: "Tien" or "Hau"
    
    Returns:
        DataFrame with time_s, Q_m3s columns
    """
    n_steps = int(duration_days * 86400 / DT_DISCHARGE) + 1
    time_s = np.arange(n_steps) * DT_DISCHARGE
    
    # Get fraction
    fraction = TIEN_FRACTION if river == "Tien" else HAU_FRACTION
    
    Q_base = Q_TOTAL_MEAN * fraction
    
    # Add daily variation (peak at noon, low at midnight)
    daily_cycle = Q_TOTAL_AMPLITUDE * fraction * np.cos(2 * np.pi * time_s / 86400 - np.pi)
    
    # Add some noise (±3%)
    np.random.seed(42 if river == "Tien" else 43)
    noise = Q_base * 0.03 * np.random.randn(n_steps)
    
    Q = Q_base + daily_cycle + noise
    Q = np.maximum(Q, Q_base * 0.7)  # Ensure positive
    
    return pd.DataFrame({"time_s": time_s, "Q_m3s": Q})


def generate_tide(duration_days: float, mouth: str) -> pd.DataFrame:
    """
    Generate tidal elevation time series.
    
    Args:
        duration_days: Simulation duration [days]
        mouth: One of "MyTho", "HamLuong", "CoChien", "DinhAn", "TranDe"
    
    Returns:
        DataFrame with time_s, H_m columns
    """
    params = TIDAL_MOUTHS[mouth]
    amplitude = params["amplitude"]
    phase_rad = np.deg2rad(params["phase_deg"])
    
    n_steps = int(duration_days * 86400 / DT_TIDE) + 1
    time_s = np.arange(n_steps) * DT_TIDE
    
    # Spring-neap modulation
    spring_neap = 1.0 + MSF_AMPLITUDE * np.cos(MSF_OMEGA * time_s)
    
    # M2 tide with phase
    H = amplitude * spring_neap * np.cos(M2_OMEGA * time_s - phase_rad)
    
    # Add small noise (±2 cm)
    np.random.seed(hash(mouth) % 2**31)
    noise = 0.02 * np.random.randn(n_steps)
    H = H + noise
    
    return pd.DataFrame({"time_s": time_s, "H_m": H})


def generate_species_timeseries(duration_days: float, 
                                 species_dict: dict,
                                 noise_pct: float = 0.05) -> pd.DataFrame:
    """
    Generate time-varying species concentrations.
    
    Args:
        duration_days: Simulation duration [days]
        species_dict: Dictionary of species name -> concentration
        noise_pct: Random noise percentage
    
    Returns:
        DataFrame with time_s column and species columns
    """
    n_steps = int(duration_days * 86400 / DT_SPECIES) + 1
    time_s = np.arange(n_steps) * DT_SPECIES
    
    data = {"time_s": time_s}
    
    np.random.seed(12345)
    for name, conc in species_dict.items():
        # Add small random variation
        noise = 1.0 + noise_pct * np.random.randn(n_steps)
        noise = np.clip(noise, 0.9, 1.1)
        data[name] = conc * noise
    
    return pd.DataFrame(data)


# ===========================================================================
# MAIN EXECUTION
# ===========================================================================

def main():
    """Generate all forcing data files."""
    
    print("=" * 70)
    print("Generating Synthetic Forcing Data: Mekong Delta Full")
    print("=" * 70)
    print(f"Start date: {START_DATE}")
    print(f"Duration: {DURATION_DAYS} days + {WARMUP_DAYS} days warmup = {TOTAL_DAYS} days")
    print(f"Output directory: {FORCING_DIR}")
    print()
    
    # Create directories
    FORCING_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # ========================================================================
    # 1. DISCHARGE FILES
    # ========================================================================
    print("1. Generating discharge time series...")
    
    # Tien Inlet
    Q_tien = generate_discharge(TOTAL_DAYS, "Tien")
    Q_tien.to_csv(FORCING_DIR / "Tien_Inlet.csv", index=False)
    print(f"   - Tien_Inlet.csv: Q_mean = {Q_tien['Q_m3s'].mean():.0f} m³/s")
    
    # Hau Inlet
    Q_hau = generate_discharge(TOTAL_DAYS, "Hau")
    Q_hau.to_csv(FORCING_DIR / "Hau_Inlet.csv", index=False)
    print(f"   - Hau_Inlet.csv: Q_mean = {Q_hau['Q_m3s'].mean():.0f} m³/s")
    
    # ========================================================================
    # 2. TIDAL BOUNDARY FILES
    # ========================================================================
    print("\n2. Generating tidal boundary conditions...")
    
    for mouth in ["MyTho", "HamLuong", "CoChien", "DinhAn", "TranDe", "Hau"]:
        tide = generate_tide(TOTAL_DAYS, mouth)
        filename = f"{mouth}_Tide.csv"
        tide.to_csv(FORCING_DIR / filename, index=False)
        params = TIDAL_MOUTHS[mouth]
        print(f"   - {filename}: amp = {params['amplitude']:.1f} m, "
              f"phase = {params['phase_deg']:.0f}°")
    
    # ========================================================================
    # 3. SPECIES BOUNDARY FILES (with time_s column!)
    # ========================================================================
    print("\n3. Generating species boundary conditions...")
    
    # River species (upstream)
    species_river = generate_species_timeseries(TOTAL_DAYS, SPECIES_RIVER)
    species_river.to_csv(FORCING_DIR / "species_river_ts.csv", index=False)
    print(f"   - species_river_ts.csv: {len(SPECIES_RIVER)} species")
    print(f"     TOC = {SPECIES_RIVER['toc']:.0f} µmolC/L, O2 = {SPECIES_RIVER['o2']:.0f} µmol/L")
    
    # Ocean species (downstream)
    species_ocean = generate_species_timeseries(TOTAL_DAYS, SPECIES_OCEAN)
    species_ocean.to_csv(FORCING_DIR / "species_ocean_ts.csv", index=False)
    print(f"   - species_ocean_ts.csv: {len(SPECIES_OCEAN)} species")
    print(f"     Salinity = {SPECIES_OCEAN['salinity']:.0f} PSU, O2 = {SPECIES_OCEAN['o2']:.0f} µmol/L")
    
    # ========================================================================
    # 4. SUMMARY
    # ========================================================================
    print("\n" + "=" * 70)
    print("Forcing Data Generation Complete!")
    print("=" * 70)
    print()
    print("Files created:")
    print("  Discharge:")
    print(f"    - {FORCING_DIR}/Tien_Inlet.csv")
    print(f"    - {FORCING_DIR}/Hau_Inlet.csv")
    print("  Tides:")
    for mouth in TIDAL_MOUTHS:
        print(f"    - {FORCING_DIR}/{mouth}_Tide.csv")
    print("  Species (time-series):")
    print(f"    - {FORCING_DIR}/species_river_ts.csv")
    print(f"    - {FORCING_DIR}/species_ocean_ts.csv")
    print()
    print("Hydrological summary:")
    print(f"  Total Mekong Q at border: ~{Q_TOTAL_MEAN:.0f} m³/s (dry season)")
    print(f"  Tien River: {TIEN_FRACTION*100:.0f}% = ~{Q_TOTAL_MEAN*TIEN_FRACTION:.0f} m³/s")
    print(f"  Hau River: {HAU_FRACTION*100:.0f}% = ~{Q_TOTAL_MEAN*HAU_FRACTION:.0f} m³/s")
    print()
    print("Next step: Run run_validation.py to test the model")
    
    return 0


if __name__ == "__main__":
    exit(main())
