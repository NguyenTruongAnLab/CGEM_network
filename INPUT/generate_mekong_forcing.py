#!/usr/bin/env python3
"""
Mekong Delta Forcing Generator for Parallel Estuary Models
Scientific approach: Savenije (2012) + Jordan et al. (2019) flow split

Generates:
  - Hau_Input.csv: 40% of total discharge
  - Tien_Input.csv: 60% of total discharge
  - Tidal boundary conditions for all mouths
  
Flow characteristics:
  - Dry season base: 2,000 m³/s
  - Flood peak (September): 35,000 m³/s
  - Annual mean: 14,600 m³/s
  
Citation: Jordan et al. (2019) Nature Scientific Reports
          MRC (Mekong River Commission) gauge data
"""

from datetime import datetime, timedelta
import os
import math
import random

# ============================================================================
# CONFIGURATION
# ============================================================================

YEAR = 2017
OUTPUT_DIR = "forcing_data"

# Discharge parameters (m³/s) - Based on MRC official data
Q_BASE = 2000.0          # Dry season base flow
Q_FLOOD_PEAK = 33000.0   # Peak flood amplitude (added to base)

# Flow split - Based on Jordan et al. (2019) and MRC observations
# "Tien carries ~60% of annual flow, Hau ~40%"
RATIO_TIEN = 0.60
RATIO_HAU = 0.40

# Tidal parameters - M2 (semi-diurnal) + K1 (diurnal) components
# Amplitudes from Nguyen et al. (2006) field measurements
TIDAL_STATIONS = {
    # Hau Combined: Aggregate of Dinh An + Tran De mouths
    'Hau_Combined_Tide.csv': {
        'mean': 0.0, 
        'A_M2': 2.1,   # Semi-diurnal amplitude (m)
        'A_K1': 0.7,   # Diurnal amplitude (m)
        'phi_M2': 0.0, 
        'phi_K1': 0.0
    },
    
    # Tien distributaries - Individual mouth characteristics
    'HamLuong_Tide.csv': {
        'mean': 0.0,
        'A_M2': 1.7,
        'A_K1': 0.5,
        'phi_M2': -0.3,
        'phi_K1': 0.2
    },
    
    'CoChien_Tide.csv': {
        'mean': 0.0,
        'A_M2': 2.1,
        'A_K1': 0.7,
        'phi_M2': 0.2,
        'phi_K1': -0.1
    },
    
    'MyTho_Tide.csv': {
        'mean': 0.0,
        'A_M2': 1.6,
        'A_K1': 0.5,
        'phi_M2': 0.5,
        'phi_K1': 0.3
    }
}

# Tidal periods (hours)
M2_PERIOD = 12.42  # Principal lunar semi-diurnal
K1_PERIOD = 23.93  # Lunar diurnal
SPRING_NEAP_PERIOD = 14.77 * 24  # Spring-neap cycle (hours)

# Random seed for reproducibility
RANDOM_SEED = 42

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def ensure_directory(path):
    """Create directory if it doesn't exist"""
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Created directory: {path}")

def flood_pulse_shape(day_of_year):
    """
    Generate Gaussian flood pulse centered on day 265 (late September)
    
    Peak monsoon in Mekong occurs August-September
    This creates a realistic annual hydrograph
    
    Returns: Multiplier [0, 1] representing flood intensity
    """
    # Gaussian centered at day 265 with sigma = 70 days
    # This gives ~4 month flood season (July-October)
    return math.exp(-((day_of_year - 265)**2) / (70**2))

# ============================================================================
# DISCHARGE GENERATION
# ============================================================================

def write_discharge_file(filename, split_ratio):
    """
    Generate daily discharge time series for one year
    
    Args:
        filename: Output CSV filename
        split_ratio: Fraction of total discharge for this branch
    """
    filepath = os.path.join(OUTPUT_DIR, filename)
    print(f"Generating {filename} (Flow split: {split_ratio*100:.0f}%)...")
    
    random.seed(RANDOM_SEED)
    
    with open(filepath, 'w') as f:
        # Write header (required by model)
        f.write('Date,Q\\n')
        
        # Generate daily values for entire year
        current_date = datetime(YEAR, 1, 1)
        end_date = datetime(YEAR, 12, 31)
        
        while current_date <= end_date:
            doy = current_date.timetuple().tm_yday
            
            # Calculate total discharge with seasonal variation
            pulse = flood_pulse_shape(doy)
            Q_total = Q_BASE + Q_FLOOD_PEAK * pulse
            
            # Add hydrologic variability (±5% noise)
            noise = random.gauss(0, Q_total * 0.05)
            Q_total_with_noise = Q_total + noise
            
            # Ensure non-negative (physical constraint)
            Q_total_with_noise = max(Q_total_with_noise, 500.0)
            
            # Apply flow split for this branch
            Q_branch = Q_total_with_noise * split_ratio
            
            # Write to file
            f.write(f"{current_date.strftime('%Y-%m-%d')},{Q_branch:.2f}\\n")
            
            current_date += timedelta(days=1)
    
    print(f"  ✓ Generated {filename} (365 days)")

# ============================================================================
# TIDAL GENERATION
# ============================================================================

def write_tidal_file(filename, params):
    """
    Generate hourly tidal elevation time series for one year
    
    Uses harmonic analysis: H(t) = mean + M2*sin() + K1*sin() + spring-neap
    
    Args:
        filename: Output CSV filename
        params: Dictionary with tidal constituent amplitudes and phases
    """
    filepath = os.path.join(OUTPUT_DIR, filename)
    print(f"Generating {filename}...")
    
    random.seed(RANDOM_SEED)
    
    with open(filepath, 'w') as f:
        # Write header (required by model)
        f.write('DateTime,H\\n')
        
        # Generate hourly values for entire year
        current_time = datetime(YEAR, 1, 1, 0, 0)
        end_time = datetime(YEAR, 12, 31, 23, 0)
        hour_index = 0
        
        while current_time <= end_time:
            # Calculate tidal elevation from harmonic constituents
            h = params['mean']
            
            # M2 component (semi-diurnal, ~12.42 hour period)
            h += params['A_M2'] * math.sin(
                2 * math.pi * hour_index / M2_PERIOD + params['phi_M2']
            )
            
            # K1 component (diurnal, ~24 hour period)
            h += params['A_K1'] * math.sin(
                2 * math.pi * hour_index / K1_PERIOD + params['phi_K1']
            )
            
            # Spring-neap modulation (~14.77 day cycle)
            # Amplitude varies ±10% on fortnightly cycle
            spring_neap_factor = 1.0 + 0.1 * math.sin(
                2 * math.pi * hour_index / SPRING_NEAP_PERIOD
            )
            h *= spring_neap_factor
            
            # Add small atmospheric/weather noise (±5 cm)
            h += random.gauss(0, 0.05)
            
            # Write to file
            f.write(f"{current_time.strftime('%Y-%m-%d %H:%M')},{h:.4f}\\n")
            
            current_time += timedelta(hours=1)
            hour_index += 1
    
    print(f"  ✓ Generated {filename} (8760 hours)")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """
    Main execution: Generate all forcing files for both cases
    """
    print("=" * 70)
    print("MEKONG DELTA FORCING GENERATOR")
    print("Parallel Estuary Approach - Savenije (2012)")
    print("=" * 70)
    print()
    
    # Create output directory
    ensure_directory(OUTPUT_DIR)
    print()
    
    # Generate discharge files
    print("STEP 1: Generating discharge forcing files")
    print("-" * 70)
    write_discharge_file('Hau_Input.csv', RATIO_HAU)
    write_discharge_file('Tien_Input.csv', RATIO_TIEN)
    print()
    
    # Generate tidal files
    print("STEP 2: Generating tidal forcing files")
    print("-" * 70)
    for filename, params in TIDAL_STATIONS.items():
        write_tidal_file(filename, params)
    print()
    
    # Summary
    print("=" * 70)
    print("GENERATION COMPLETE")
    print("=" * 70)
    print(f"Output directory: {OUTPUT_DIR}/")
    print()
    print("Files generated:")
    print("  Discharge (daily):")
    print("    - Hau_Input.csv (40% of total)")
    print("    - Tien_Input.csv (60% of total)")
    print()
    print("  Tidal elevation (hourly):")
    print("    - Hau_Combined_Tide.csv")
    print("    - HamLuong_Tide.csv")
    print("    - CoChien_Tide.csv")
    print("    - MyTho_Tide.csv")
    print()
    print("Ready for simulation:")
    print("  - Case Hau_River")
    print("  - Case Tien_River")
    print("=" * 70)

if __name__ == '__main__':
    main()