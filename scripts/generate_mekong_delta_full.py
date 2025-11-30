#!/usr/bin/env python3
"""
Generate Forcing Data for Full Mekong Delta Case
=================================================

This script generates all boundary condition forcing files for the 
complete Mekong Delta network (Tien + Hau + Vam Nao).

Outputs:
- forcing_data/Tien_Inlet.csv   (80% of total Mekong discharge)
- forcing_data/Hau_Inlet.csv    (20% of total Mekong discharge)
- forcing_data/MyTho_Tide.csv   (phase = 0°)
- forcing_data/HamLuong_Tide.csv (phase = 10°)
- forcing_data/CoChien_Tide.csv  (phase = 20°)
- forcing_data/DinhAn_Tide.csv   (phase = 30°)
- forcing_data/TranDe_Tide.csv   (phase = 40°)

Usage:
    python generate_mekong_delta_full.py
    
    Or use generate_mekong_inputs.py with --case Mekong_Delta_Full

Author: CGEM Development Team
Date: November 2024
"""

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from pathlib import Path

# ============================================================================
# MEKONG DELTA HYDROLOGY PARAMETERS
# ============================================================================

# Total Mekong discharge at Vietnamese border (m³/s)
# Based on MRC data at Tan Chau + Chau Doc combined
MEKONG_DISCHARGE = {
    'dry_season': {         # Dec-Apr
        'min': 4000,
        'max': 7000,
        'mean': 5500
    },
    'wet_season': {         # Jul-Oct  
        'min': 15000,
        'max': 35000,
        'mean': 22000
    },
    'transition': {         # May-Jun, Nov
        'min': 7000,
        'max': 15000,
        'mean': 10000
    }
}

# Tien/Hau split at border (before Vam Nao)
# Historical: ~80% Tien, ~20% Hau at border
# Vam Nao redistributes to achieve ~55% Tien, ~45% Hau at delta
DISCHARGE_SPLIT = {
    'Tien_fraction': 0.80,  # 80% enters via Tien at Tan Chau
    'Hau_fraction': 0.20    # 20% enters via Hau at Chau Doc
}

# Tidal parameters for each mouth
# M2 dominant with phase shifts (tide propagates SW)
TIDAL_MOUTHS = {
    'MyTho': {
        'amplitude': 1.8,   # M2 amplitude [m]
        'phase_deg': 0,     # Reference phase
        'damping': 1.0      # Damping factor
    },
    'HamLuong': {
        'amplitude': 1.9,
        'phase_deg': 8,     # 8 degrees behind My Tho
        'damping': 1.0
    },
    'CoChien': {
        'amplitude': 2.0,
        'phase_deg': 15,    # 15 degrees behind
        'damping': 1.0
    },
    'DinhAn': {
        'amplitude': 2.2,   # Higher amplitude on Hau
        'phase_deg': 25,
        'damping': 1.0
    },
    'TranDe': {
        'amplitude': 2.0,
        'phase_deg': 35,    # Most SW, largest phase lag
        'damping': 1.0
    }
}


def get_season(date: datetime) -> str:
    """Determine Mekong season from date."""
    month = date.month
    if month in [12, 1, 2, 3, 4]:
        return 'dry_season'
    elif month in [7, 8, 9, 10]:
        return 'wet_season'
    else:
        return 'transition'


def seasonal_discharge_factor(date: datetime) -> float:
    """Calculate smooth seasonal factor for discharge (0-1)."""
    doy = date.timetuple().tm_yday
    # Wet season peaks around day 245 (early September)
    wet_phase = 2 * np.pi * (doy - 245) / 365
    return 0.5 * (1 + np.cos(wet_phase))


def generate_discharge(start_date: str, duration_days: int, 
                       river: str = 'Tien') -> pd.DataFrame:
    """
    Generate discharge time series for Tien or Hau inlet.
    
    Args:
        start_date: Start date (YYYY-MM-DD)
        duration_days: Duration in days
        river: 'Tien' or 'Hau'
    
    Returns:
        DataFrame with time_s, Q_m3s columns
    """
    start = datetime.strptime(start_date, '%Y-%m-%d')
    dt_hours = 1.0  # Hourly discharge
    n_steps = int(duration_days * 24 / dt_hours) + 1
    
    # Get split fraction
    if river == 'Tien':
        fraction = DISCHARGE_SPLIT['Tien_fraction']
    else:
        fraction = DISCHARGE_SPLIT['Hau_fraction']
    
    time_s = []
    Q_m3s = []
    
    for i in range(n_steps):
        t = start + timedelta(hours=i * dt_hours)
        season = get_season(t)
        sf = seasonal_discharge_factor(t)
        
        # Base discharge from seasonal range
        Q_range = MEKONG_DISCHARGE[season]
        Q_seasonal = Q_range['min'] + sf * (Q_range['max'] - Q_range['min'])
        
        # Apply river split
        Q_river = Q_seasonal * fraction
        
        # Add daily variation (~5%)
        daily_var = 1.0 + 0.05 * np.sin(2 * np.pi * t.hour / 24)
        
        # Add random noise (~3%)
        noise = 1.0 + 0.03 * np.random.randn()
        
        Q_final = Q_river * daily_var * max(0.8, min(1.2, noise))
        
        time_s.append(i * dt_hours * 3600)
        Q_m3s.append(Q_final)
    
    return pd.DataFrame({'time_s': time_s, 'Q_m3s': Q_m3s})


def generate_tide(start_date: str, duration_days: int,
                  mouth: str) -> pd.DataFrame:
    """
    Generate M2 tidal elevation time series for a river mouth.
    
    Args:
        start_date: Start date (YYYY-MM-DD)
        duration_days: Duration in days
        mouth: One of 'MyTho', 'HamLuong', 'CoChien', 'DinhAn', 'TranDe'
    
    Returns:
        DataFrame with time_s, H_m columns
    """
    dt_hours = 0.5  # Half-hourly for tide resolution
    n_steps = int(duration_days * 24 / dt_hours) + 1
    
    # Get mouth parameters
    params = TIDAL_MOUTHS[mouth]
    amplitude = params['amplitude']
    phase_deg = params['phase_deg']
    damping = params['damping']
    
    # M2 tidal parameters
    M2_period_s = 12.42 * 3600  # M2 period in seconds
    M2_omega = 2 * np.pi / M2_period_s
    
    # Spring-neap modulation (Msf)
    Msf_period_s = 14.77 * 86400
    Msf_omega = 2 * np.pi / Msf_period_s
    
    # Phase in radians
    phase_rad = np.deg2rad(phase_deg)
    
    time_s = []
    H_m = []
    
    for i in range(n_steps):
        ts = i * dt_hours * 3600
        
        # Spring-neap modulation (±20%)
        spring_neap = 1.0 + 0.2 * np.cos(Msf_omega * ts)
        
        # M2 tide with phase shift
        H_tide = amplitude * damping * spring_neap * np.cos(M2_omega * ts - phase_rad)
        
        # Add small noise (~2 cm)
        noise = np.random.randn() * 0.02
        
        time_s.append(ts)
        H_m.append(H_tide + noise)
    
    return pd.DataFrame({'time_s': time_s, 'H_m': H_m})


def main():
    """Generate all forcing files for Full Mekong Delta case."""
    
    # Configuration
    start_date = '2024-03-01'  # Dry season start
    duration_days = 45         # Match case_config.txt (30 days + 15 warmup)
    
    # Output directory
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    forcing_dir = project_root / 'forcing_data'
    forcing_dir.mkdir(exist_ok=True)
    
    print("="*70)
    print("Generating Forcing Data: Full Mekong Delta Case")
    print("="*70)
    print(f"Start: {start_date}")
    print(f"Duration: {duration_days} days")
    print(f"Output: {forcing_dir}")
    print()
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # ========================================================================
    # 1. DISCHARGE FILES
    # ========================================================================
    print("1. Generating discharge time series...")
    
    # Tien Inlet (Tan Chau) - 80% of total
    Q_tien = generate_discharge(start_date, duration_days, river='Tien')
    Q_tien.to_csv(forcing_dir / 'Tien_Inlet.csv', index=False)
    print(f"   - Tien_Inlet.csv: Q_mean={Q_tien['Q_m3s'].mean():.0f} m³/s")
    
    # Hau Inlet (Chau Doc) - 20% of total
    Q_hau = generate_discharge(start_date, duration_days, river='Hau')
    Q_hau.to_csv(forcing_dir / 'Hau_Inlet.csv', index=False)
    print(f"   - Hau_Inlet.csv: Q_mean={Q_hau['Q_m3s'].mean():.0f} m³/s")
    
    # ========================================================================
    # 2. TIDAL BOUNDARY FILES (5 mouths)
    # ========================================================================
    print("\n2. Generating tidal boundary conditions...")
    
    for mouth in ['MyTho', 'HamLuong', 'CoChien', 'DinhAn', 'TranDe']:
        tide = generate_tide(start_date, duration_days, mouth)
        filename = f'{mouth}_Tide.csv'
        tide.to_csv(forcing_dir / filename, index=False)
        print(f"   - {filename}: amp={TIDAL_MOUTHS[mouth]['amplitude']:.1f}m, "
              f"phase={TIDAL_MOUTHS[mouth]['phase_deg']}°")
    
    # ========================================================================
    # 3. SUMMARY
    # ========================================================================
    print("\n" + "="*70)
    print("Forcing data generation complete!")
    print("="*70)
    print("\nFiles created:")
    print("  Discharge:")
    print("    - forcing_data/Tien_Inlet.csv")
    print("    - forcing_data/Hau_Inlet.csv")
    print("  Tides:")
    print("    - forcing_data/MyTho_Tide.csv")
    print("    - forcing_data/HamLuong_Tide.csv")
    print("    - forcing_data/CoChien_Tide.csv")
    print("    - forcing_data/DinhAn_Tide.csv")
    print("    - forcing_data/TranDe_Tide.csv")
    print("\nNote: Species BCs use static files in case folder:")
    print("  - INPUT/Cases/Mekong_Delta_Full/species_river.csv")
    print("  - INPUT/Cases/Mekong_Delta_Full/species_ocean.csv")
    print("\nRun simulation with:")
    print("  ./bin/Debug/CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt")


if __name__ == '__main__':
    main()
