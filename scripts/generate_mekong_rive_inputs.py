#!/usr/bin/env python3
"""
CGEM-RIVE Mekong Delta Input Data Generator

Generates realistic biogeochemical input data for the Mekong Delta based on:
- MRC (Mekong River Commission) monitoring data
- Published literature values
- Seasonal patterns typical of tropical monsoon estuaries

Output files are CSV format compatible with CGEM_network.

Reference ranges based on:
- Nguyen et al. (2019) Mekong biogeochemistry
- Toming et al. (2020) DOC in tropical rivers
- Volta et al. (2016) Scheldt estuary pCO2
- Billen et al. (1994) RIVE model parameters

Author: CGEM-RIVE Development Team
Date: 2024
"""

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import os
import argparse

# ============================================================================
# MEKONG DELTA TYPICAL RANGES (from literature)
# ============================================================================

MEKONG_RANGES = {
    # Physical parameters
    'water_temp': (25.0, 32.0),          # [°C] - tropical, seasonal variation
    'salinity_ocean': (28.0, 35.0),      # [PSU] - South China Sea
    'salinity_river': (0.0, 0.5),        # [PSU] - freshwater from upstream
    
    # Nutrients (µmol/L)
    'no3_river': (20.0, 80.0),           # Higher in wet season (agriculture runoff)
    'no3_ocean': (0.5, 5.0),             # Low in coastal ocean
    'nh4_river': (5.0, 30.0),            # Agricultural/urban sources
    'nh4_ocean': (0.5, 2.0),             # Low in ocean
    'po4_river': (0.5, 3.0),             # Fertilizer runoff
    'po4_ocean': (0.1, 0.5),             # Low in ocean
    'dsi_river': (100.0, 200.0),         # High silica from weathering
    'dsi_ocean': (5.0, 20.0),            # Low in ocean
    
    # Organic carbon (µmol C/L)
    'doc_river': (200.0, 500.0),         # High DOC in tropical rivers
    'doc_ocean': (50.0, 100.0),          # Lower in ocean
    'poc_river': (100.0, 300.0),         # Variable with SPM
    'poc_ocean': (20.0, 50.0),           # Lower in ocean
    
    # Dissolved oxygen (µmol/L)
    'o2_river': (200.0, 280.0),          # Near saturation
    'o2_ocean': (200.0, 250.0),          # Well-oxygenated
    
    # Carbonate system
    'dic_river': (1800.0, 2400.0),       # [µmol/L] - high from weathering
    'dic_ocean': (2000.0, 2200.0),       # [µmol/L] - oceanic
    'at_river': (1800.0, 2600.0),        # [µmol/L] - alkalinity
    'at_ocean': (2200.0, 2400.0),        # [µmol/L] - oceanic
    
    # Suspended particulate matter
    'spm_river': (50.0, 300.0),          # [mg/L] - high in wet season
    'spm_ocean': (5.0, 30.0),            # [mg/L] - lower in ocean
    
    # Phytoplankton (mg C/L)
    'phy_river': (0.1, 1.0),             # Low upstream (turbidity)
    'phy_ocean': (0.5, 3.0),             # Higher at coast (light)
    
    # RIVE-specific: Bacteria (mg C/L)
    'bag_init': (0.05, 0.2),             # Attached bacteria
    'bap_init': (0.02, 0.1),             # Free bacteria
    
    # RIVE-specific: DSS (mg C/L)
    'dss_init': (0.5, 2.0),              # Dissolved simple substrates
    
    # Discharge (m³/s) - Tien River
    'Q_dry': (1500.0, 3000.0),           # Dec-Apr
    'Q_wet': (8000.0, 15000.0),          # Jul-Oct
    'Q_transition': (3000.0, 8000.0),    # May-Jun, Nov
    
    # Tidal amplitude (m)
    'tide_amplitude': (1.5, 2.5),        # Spring-neap cycle
}


def seasonal_factor(day_of_year: int) -> dict:
    """
    Calculate seasonal factors based on day of year.
    Mekong monsoon: wet season Jul-Oct, dry season Dec-Apr
    
    Returns dict with factors for different parameters.
    """
    # Wet season peaks around day 240 (late August)
    wet_phase = 2 * np.pi * (day_of_year - 240) / 365
    wet_factor = 0.5 * (1 + np.cos(wet_phase))  # 1 in wet, 0 in dry
    
    # Temperature peaks in April (day 100)
    temp_phase = 2 * np.pi * (day_of_year - 100) / 365
    temp_factor = 0.5 * (1 + np.cos(temp_phase))
    
    return {
        'discharge': 0.3 + 0.7 * wet_factor,      # Higher in wet season
        'spm': 0.2 + 0.8 * wet_factor,            # Higher turbidity in wet
        'nutrients': 0.4 + 0.6 * wet_factor,      # More runoff in wet
        'temperature': 0.6 + 0.4 * temp_factor,   # Warmer in Apr-May
        'salinity': 1.0 - 0.6 * wet_factor,       # Lower intrusion in wet
    }


def interpolate_range(low: float, high: float, factor: float) -> float:
    """Interpolate between low and high based on factor [0-1]."""
    return low + factor * (high - low)


def generate_discharge_timeseries(
    start_date: str,
    duration_days: int,
    dt_hours: float = 1.0,
    river: str = 'Tien'
) -> pd.DataFrame:
    """
    Generate realistic Mekong discharge time series.
    
    Args:
        start_date: Start date (YYYY-MM-DD)
        duration_days: Duration in days
        dt_hours: Time step in hours
        river: 'Tien' or 'Hau'
        
    Returns:
        DataFrame with columns: time_s, Q_m3s
    """
    start = datetime.strptime(start_date, '%Y-%m-%d')
    n_steps = int(duration_days * 24 / dt_hours) + 1
    
    times = [start + timedelta(hours=i * dt_hours) for i in range(n_steps)]
    time_seconds = np.arange(n_steps) * dt_hours * 3600
    
    # Scale factor for Tien vs Hau (Tien is ~60% of total)
    scale = 1.0 if river == 'Tien' else 0.67
    
    Q_values = []
    for t in times:
        doy = t.timetuple().tm_yday
        sf = seasonal_factor(doy)
        
        # Base discharge with seasonal variation
        Q_dry = np.mean(MEKONG_RANGES['Q_dry'])
        Q_wet = np.mean(MEKONG_RANGES['Q_wet'])
        Q_base = Q_dry + sf['discharge'] * (Q_wet - Q_dry)
        
        # Add daily variation (diurnal discharge pattern from upstream)
        hour_factor = 1.0 + 0.05 * np.sin(2 * np.pi * t.hour / 24)
        
        # Add noise (weather variability)
        noise = 1.0 + 0.1 * np.random.randn()
        
        Q = Q_base * scale * hour_factor * max(0.5, noise)
        Q_values.append(Q)
    
    return pd.DataFrame({
        'time_s': time_seconds,
        'Q_m3s': Q_values
    })


def generate_tide_timeseries(
    start_date: str,
    duration_days: int,
    dt_hours: float = 0.5,
    location: str = 'mouth'
) -> pd.DataFrame:
    """
    Generate M2 tidal elevation time series.
    
    Args:
        start_date: Start date (YYYY-MM-DD)
        duration_days: Duration in days
        dt_hours: Time step in hours
        location: 'mouth', 'mid', or 'upstream'
        
    Returns:
        DataFrame with columns: time_s, H_m
    """
    n_steps = int(duration_days * 24 / dt_hours) + 1
    time_seconds = np.arange(n_steps) * dt_hours * 3600
    
    # M2 tidal period (12.42 hours)
    M2_period = 12.42 * 3600  # seconds
    M2_freq = 2 * np.pi / M2_period
    
    # Spring-neap modulation (14.77 days)
    Msf_period = 14.77 * 86400  # seconds
    Msf_freq = 2 * np.pi / Msf_period
    
    # Amplitude depends on location (dampens upstream)
    amp_base = np.mean(MEKONG_RANGES['tide_amplitude'])
    if location == 'mouth':
        amp = amp_base
    elif location == 'mid':
        amp = amp_base * 0.7
    else:  # upstream
        amp = amp_base * 0.4
    
    # Mean sea level reference
    msl = 0.0
    
    H_values = []
    for ts in time_seconds:
        # Spring-neap modulation (±20%)
        spring_neap = 1.0 + 0.2 * np.cos(Msf_freq * ts)
        
        # M2 tide
        H = msl + amp * spring_neap * np.cos(M2_freq * ts)
        H_values.append(H)
    
    return pd.DataFrame({
        'time_s': time_seconds,
        'H_m': H_values
    })


def generate_species_boundary(
    start_date: str,
    duration_days: int,
    dt_hours: float = 6.0,
    boundary_type: str = 'river'
) -> pd.DataFrame:
    """
    Generate biogeochemical boundary condition time series.
    
    Args:
        start_date: Start date (YYYY-MM-DD)
        duration_days: Duration in days
        dt_hours: Time step in hours
        boundary_type: 'river' or 'ocean'
        
    Returns:
        DataFrame with all species concentrations
    """
    start = datetime.strptime(start_date, '%Y-%m-%d')
    n_steps = int(duration_days * 24 / dt_hours) + 1
    
    times = [start + timedelta(hours=i * dt_hours) for i in range(n_steps)]
    time_seconds = np.arange(n_steps) * dt_hours * 3600
    
    is_river = (boundary_type == 'river')
    
    data = {'time_s': time_seconds}
    
    for t_idx, t in enumerate(times):
        doy = t.timetuple().tm_yday
        sf = seasonal_factor(doy)
        
        # Add noise
        def noisy(val, noise_pct=0.1):
            return val * (1.0 + noise_pct * np.random.randn())
        
        # Salinity
        if is_river:
            sal = interpolate_range(*MEKONG_RANGES['salinity_river'], np.random.rand() * 0.3)
        else:
            sal = interpolate_range(*MEKONG_RANGES['salinity_ocean'], 0.5 + 0.5 * sf['salinity'])
        data.setdefault('salinity', []).append(noisy(sal, 0.05))
        
        # Temperature
        temp_range = MEKONG_RANGES['water_temp']
        temp = interpolate_range(*temp_range, sf['temperature'])
        data.setdefault('temperature', []).append(noisy(temp, 0.02))
        
        # Nutrients
        no3_range = MEKONG_RANGES['no3_river'] if is_river else MEKONG_RANGES['no3_ocean']
        nh4_range = MEKONG_RANGES['nh4_river'] if is_river else MEKONG_RANGES['nh4_ocean']
        po4_range = MEKONG_RANGES['po4_river'] if is_river else MEKONG_RANGES['po4_ocean']
        dsi_range = MEKONG_RANGES['dsi_river'] if is_river else MEKONG_RANGES['dsi_ocean']
        
        data.setdefault('no3', []).append(noisy(interpolate_range(*no3_range, sf['nutrients'])))
        data.setdefault('nh4', []).append(noisy(interpolate_range(*nh4_range, sf['nutrients'])))
        data.setdefault('po4', []).append(noisy(interpolate_range(*po4_range, sf['nutrients'])))
        data.setdefault('dsi', []).append(noisy(interpolate_range(*dsi_range, 0.5)))
        
        # Organic carbon
        doc_range = MEKONG_RANGES['doc_river'] if is_river else MEKONG_RANGES['doc_ocean']
        poc_range = MEKONG_RANGES['poc_river'] if is_river else MEKONG_RANGES['poc_ocean']
        
        doc = noisy(interpolate_range(*doc_range, sf['nutrients']))
        poc = noisy(interpolate_range(*poc_range, sf['spm']))
        
        # Partition DOC into RIVE pools
        frac_hd1 = 0.15  # Labile
        frac_hd2 = 0.35  # Semi-labile
        # Rest is HD3 (refractory)
        
        data.setdefault('hd1', []).append(doc * frac_hd1)
        data.setdefault('hd2', []).append(doc * frac_hd2)
        data.setdefault('hd3', []).append(doc * (1.0 - frac_hd1 - frac_hd2))
        
        # Partition POC into RIVE pools
        frac_hp1 = 0.10
        frac_hp2 = 0.30
        
        data.setdefault('hp1', []).append(poc * frac_hp1)
        data.setdefault('hp2', []).append(poc * frac_hp2)
        data.setdefault('hp3', []).append(poc * (1.0 - frac_hp1 - frac_hp2))
        
        # Total OC (for backward compatibility)
        data.setdefault('toc', []).append(doc + poc)
        
        # Dissolved oxygen (near saturation in river, slightly lower in ocean)
        o2_range = MEKONG_RANGES['o2_river'] if is_river else MEKONG_RANGES['o2_ocean']
        data.setdefault('o2', []).append(noisy(interpolate_range(*o2_range, 0.7)))
        
        # Carbonate system
        dic_range = MEKONG_RANGES['dic_river'] if is_river else MEKONG_RANGES['dic_ocean']
        at_range = MEKONG_RANGES['at_river'] if is_river else MEKONG_RANGES['at_ocean']
        
        data.setdefault('dic', []).append(noisy(interpolate_range(*dic_range, 0.5)))
        data.setdefault('at', []).append(noisy(interpolate_range(*at_range, 0.5)))
        
        # SPM
        spm_range = MEKONG_RANGES['spm_river'] if is_river else MEKONG_RANGES['spm_ocean']
        spm = noisy(interpolate_range(*spm_range, sf['spm']))
        data.setdefault('spm', []).append(spm)
        
        # Phytoplankton
        phy_range = MEKONG_RANGES['phy_river'] if is_river else MEKONG_RANGES['phy_ocean']
        phy = noisy(interpolate_range(*phy_range, 0.5))
        data.setdefault('phy1', []).append(phy * 0.7)  # Diatoms dominant
        data.setdefault('phy2', []).append(phy * 0.3)  # Non-siliceous
        
        # RIVE bacteria
        bag = noisy(interpolate_range(*MEKONG_RANGES['bag_init'], 0.5))
        bap = noisy(interpolate_range(*MEKONG_RANGES['bap_init'], 0.5))
        data.setdefault('bag', []).append(bag)
        data.setdefault('bap', []).append(bap)
        
        # RIVE DSS (dissolved simple substrates)
        dss = noisy(interpolate_range(*MEKONG_RANGES['dss_init'], 0.5))
        data.setdefault('dss', []).append(dss)
        
        # PIP (calculated from equilibrium with SPM)
        # Using simplified estimate: PIP ~ pac * SPM * PO4 / (kpads + PO4)
        po4_val = data['po4'][-1]
        kpads = 3.43
        pac = 0.37
        pip = pac * spm * po4_val / (kpads + po4_val) if (kpads + po4_val) > 0 else 0.0
        data.setdefault('pip', []).append(pip)
    
    return pd.DataFrame(data)


def generate_biogeo_params(output_path: str, region: str = 'mekong'):
    """
    Generate biogeo_params.txt with RIVE parameters for Mekong.
    
    Args:
        output_path: Path to write biogeo_params.txt
        region: 'mekong' or 'default'
    """
    params = f"""# =============================================================================
# CGEM-RIVE Biogeochemistry Parameters
# Region: {region.upper()}
# Generated by generate_mekong_rive_inputs.py
# =============================================================================

# -----------------------------------------------------------------------------
# Environmental parameters
# -----------------------------------------------------------------------------
water_temp = 28.0           # Mean water temperature [°C]
ws = 0.001                  # SPM settling velocity [m/s]

# -----------------------------------------------------------------------------
# Light parameters
# -----------------------------------------------------------------------------
I0 = 400.0                  # Surface irradiance [W/m²] - tropical
kd1 = 0.15                  # Base light attenuation [1/m]
kd2_spm = 0.02              # SPM light attenuation [L/mg/m]
kd2_phy1 = 0.005            # Phy1 self-shading
kd2_phy2 = 0.005            # Phy2 self-shading

# -----------------------------------------------------------------------------
# Phytoplankton parameters
# -----------------------------------------------------------------------------
# Phy1 = Diatoms (dominant in Mekong)
alpha1 = 0.02               # Initial slope [mgC/mgChl/W/m²]
pbmax1 = 2.5                # Max photosynthesis rate [mgC/mgChl/h]
kexc1 = 0.1                 # Excretion fraction [-]
kgrowth1 = 0.1              # Growth respiration fraction [-]
kmaint1 = 0.02              # Maintenance respiration [1/h]
kmort1 = 0.05               # Mortality rate [1/day]

# Phy2 = Non-siliceous algae
alpha2 = 0.015              # Initial slope
pbmax2 = 2.0                # Max photosynthesis
kexc2 = 0.1                 # Excretion fraction
kgrowth2 = 0.1              # Growth respiration
kmaint2 = 0.015             # Maintenance respiration
kmort2 = 0.04               # Mortality rate

# -----------------------------------------------------------------------------
# Nutrient limitation (half-saturation constants)
# -----------------------------------------------------------------------------
kdsi1 = 5.0                 # DSi for Phy1 [µmol/L]
kn1 = 2.0                   # N for Phy1 [µmol/L]
kpo41 = 0.5                 # PO4 for Phy1 [µmol/L]
kn2 = 3.0                   # N for Phy2 [µmol/L]
kpo42 = 0.7                 # PO4 for Phy2 [µmol/L]

# -----------------------------------------------------------------------------
# Original CGEM decomposition (used for backward compatibility)
# -----------------------------------------------------------------------------
kox = 0.1                   # Aerobic TOC decay [1/day]
kdenit = 0.05               # Denitrification [1/day]
knit = 0.1                  # Nitrification [1/day]
ktox = 50.0                 # TOC half-saturation [µmol/L]
ko2 = 10.0                  # O2 half-saturation for respiration [µmol/L]
ko2_nit = 5.0               # O2 half-saturation for nitrification [µmol/L]
kno3 = 5.0                  # NO3 half-saturation [µmol/L]
knh4 = 2.0                  # NH4 half-saturation [µmol/L]
kino2 = 5.0                 # O2 inhibition for denitrification [µmol/L]

# -----------------------------------------------------------------------------
# Stoichiometry (Redfield-like)
# -----------------------------------------------------------------------------
redn = 0.151                # N:C ratio [mol N / mol C]
redp = 0.01                 # P:C ratio [mol P / mol C]
redsi = 0.1                 # Si:C ratio for diatoms

# -----------------------------------------------------------------------------
# Gas exchange parameters
# -----------------------------------------------------------------------------
pco2_atm = 420.0            # Atmospheric pCO2 [µatm]
wind_speed = 3.5            # Mean wind at 10m [m/s]
wind_coeff = 0.251          # Wanninkhof coefficient
schmidt_exp = -0.5          # Schmidt number exponent
current_k_factor = 0.25     # Current contribution to k

# -----------------------------------------------------------------------------
# Benthic respiration (enhanced for shallow tropical estuary)
# -----------------------------------------------------------------------------
benthic_resp_20C = 80.0     # Benthic CO2 flux at 20°C [mmol C/m²/day]
benthic_Q10 = 2.0           # Temperature coefficient

# =============================================================================
# RIVE Multi-Pool Organic Matter Parameters
# Reference: Billen et al. (1994), Garnier et al. (2002)
# =============================================================================

# Hydrolysis rates [1/day] - conversion of pools
khydr1 = 0.75               # HD1/HP1 labile hydrolysis
khydr2 = 0.25               # HD2/HP2 semi-labile hydrolysis
khydr3 = 0.01               # HD3/HP3 refractory hydrolysis

# Fractionation of organic matter inputs
frac_hd1 = 0.15             # Fraction of DOC as labile HD1
frac_hd2 = 0.35             # Fraction of DOC as semi-labile HD2
frac_hp1 = 0.10             # Fraction of POC as labile HP1
frac_hp2 = 0.30             # Fraction of POC as semi-labile HP2

# =============================================================================
# RIVE Bacteria Parameters
# Reference: RIVE bacteria.cpp
# =============================================================================

# BAG = Attached (gros) bacteria - attached to particles
bag_bmax20 = 0.6            # Max growth rate at 20°C [1/h]
bag_kdb20 = 0.05            # Mortality rate at 20°C [1/h]
bag_topt = 22.0             # Optimal temperature [°C]
bag_dti = 12.0              # Temperature spread [°C]
bag_vs = 0.02               # Settling velocity [m/h]

# BAP = Free (petit) bacteria - free-floating
bap_bmax20 = 0.16           # Max growth rate at 20°C [1/h]
bap_kdb20 = 0.02            # Mortality rate at 20°C [1/h]
bap_topt = 20.0             # Optimal temperature [°C]
bap_dti = 17.0              # Temperature spread [°C]

# Common bacteria parameters
bac_ks = 0.1                # DSS half-saturation [mg C/L]
bac_yield = 0.25            # Growth yield [-]

# =============================================================================
# RIVE Phosphorus Adsorption Parameters
# Reference: RIVE phosphorus.cpp
# =============================================================================
kpads = 3.43                # Adsorption equilibrium constant
pac = 0.37                  # P adsorption capacity [µmol P / mg SPM]

# =============================================================================
# RIVE Benthic Parameters
# =============================================================================
zf_init = 0.001             # Initial fluid sediment depth [m]
benthic_porosity = 0.88     # Sediment porosity [-]
benthic_density = 2300000.0 # Sediment density [g/m³]
"""
    
    with open(output_path, 'w') as f:
        f.write(params)
    
    print(f"Generated: {output_path}")


def generate_case_inputs(
    case_dir: str,
    case_name: str,
    start_date: str = '2017-01-01',
    duration_days: int = 365
):
    """
    Generate complete set of input files for a case.
    
    Args:
        case_dir: Directory to write files
        case_name: Case name (e.g., 'Tien_River')
        start_date: Simulation start date
        duration_days: Duration in days
    """
    os.makedirs(case_dir, exist_ok=True)
    forcing_dir = os.path.join(case_dir, 'forcing_data')
    os.makedirs(forcing_dir, exist_ok=True)
    
    print(f"\nGenerating inputs for case: {case_name}")
    print(f"  Start: {start_date}, Duration: {duration_days} days")
    print(f"  Output: {case_dir}")
    
    # 1. Discharge time series
    Q_df = generate_discharge_timeseries(start_date, duration_days, dt_hours=1.0, river='Tien')
    Q_path = os.path.join(forcing_dir, 'Tien_Input.csv')
    Q_df.to_csv(Q_path, index=False)
    print(f"  - Discharge: {Q_path}")
    
    # 2. Tidal elevation time series for each branch mouth
    for branch_name in ['HamLuong', 'CoChien', 'MyTho']:
        tide_df = generate_tide_timeseries(start_date, duration_days, dt_hours=0.5, location='mouth')
        tide_path = os.path.join(forcing_dir, f'{branch_name}_Tide.csv')
        tide_df.to_csv(tide_path, index=False)
        print(f"  - Tide ({branch_name}): {tide_path}")
    
    # 3. Species boundary conditions (river)
    river_bc = generate_species_boundary(start_date, duration_days, dt_hours=6.0, boundary_type='river')
    river_bc_path = os.path.join(forcing_dir, 'species_river.csv')
    river_bc.to_csv(river_bc_path, index=False)
    print(f"  - River BC: {river_bc_path}")
    
    # 4. Species boundary conditions (ocean)
    ocean_bc = generate_species_boundary(start_date, duration_days, dt_hours=6.0, boundary_type='ocean')
    ocean_bc_path = os.path.join(forcing_dir, 'species_ocean.csv')
    ocean_bc.to_csv(ocean_bc_path, index=False)
    print(f"  - Ocean BC: {ocean_bc_path}")
    
    # 5. Biogeochemistry parameters
    biogeo_path = os.path.join(case_dir, 'biogeo_params.txt')
    generate_biogeo_params(biogeo_path, region='mekong')
    
    print(f"\nInput generation complete!")
    print(f"Total files: 6 CSV + 1 biogeo_params.txt")


def main():
    parser = argparse.ArgumentParser(
        description='Generate CGEM-RIVE input data for Mekong Delta'
    )
    parser.add_argument(
        '--case', '-c',
        default='Tien_River',
        help='Case name (default: Tien_River)'
    )
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Output directory (default: INPUT/Cases/<case_name>)'
    )
    parser.add_argument(
        '--start', '-s',
        default='2017-01-01',
        help='Start date YYYY-MM-DD (default: 2017-01-01)'
    )
    parser.add_argument(
        '--days', '-d',
        type=int,
        default=365,
        help='Duration in days (default: 365)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility'
    )
    
    args = parser.parse_args()
    
    # Set random seed
    np.random.seed(args.seed)
    
    # Determine output directory
    if args.output:
        case_dir = args.output
    else:
        # Find project root (where INPUT folder is)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
        case_dir = os.path.join(project_root, 'INPUT', 'Cases', args.case)
    
    generate_case_inputs(
        case_dir=case_dir,
        case_name=args.case,
        start_date=args.start,
        duration_days=args.days
    )


if __name__ == '__main__':
    main()
