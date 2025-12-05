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
# Output directories for each case
TIEN_FORCING_DIR = os.path.join("Cases", "Tien_River", "forcing_data")
HAU_FORCING_DIR = os.path.join("Cases", "Hau_River", "forcing_data")
# Also output to root forcing_data for convenience
ROOT_FORCING_DIR = os.path.join("..", "forcing_data")

# Discharge parameters (m³/s) - Based on MRC official data
Q_BASE = 2000.0          # Dry season base flow
Q_FLOOD_PEAK = 33000.0   # Peak flood amplitude (added to base)

# Flow split - Based on Jordan et al. (2019) and MRC observations
# "Tien carries ~60% of annual flow, Hau ~40%"
RATIO_TIEN = 0.60
RATIO_HAU = 0.40

# Tidal parameters - M2 (semi-diurnal) + K1 (diurnal) + S2 (solar) components
# Amplitudes from Nguyen et al. (2006) field measurements
# PHASE 5 UPGRADE (Dec 2025): Added S2 for proper spring-neap forcing
# Station names match boundary_map.csv in each case

# Tien River case tidal stations
TIEN_TIDAL_STATIONS = {
    'MyTho_Tide.csv': {
        'mean': 0.0,
        'A_M2': 1.3,   # Semi-diurnal amplitude (m) - reduced to account for S2
        'A_K1': 0.5,   # Diurnal amplitude (m)
        'A_S2': 0.35,  # Solar semi-diurnal (m) - NEW
        'phi_M2': 0.5,  # Phase offset
        'phi_K1': 0.3,
        'phi_S2': 0.4   # NEW
    },
    'HamLuong_Tide.csv': {
        'mean': 0.0,
        'A_M2': 1.4,
        'A_K1': 0.5,
        'A_S2': 0.38,
        'phi_M2': -0.3,
        'phi_K1': 0.2,
        'phi_S2': -0.2
    },
    'CoChien_Tide.csv': {
        'mean': 0.0,
        'A_M2': 1.7,
        'A_K1': 0.7,
        'A_S2': 0.45,
        'phi_M2': 0.2,
        'phi_K1': -0.1,
        'phi_S2': 0.1
    }
}

# Hau River case tidal station
HAU_TIDAL_STATIONS = {
    'Hau_Combined_Tide.csv': {
        'mean': 0.0, 
        'A_M2': 1.7,
        'A_K1': 0.7,
        'A_S2': 0.45,
        'phi_M2': 0.0, 
        'phi_K1': 0.0,
        'phi_S2': 0.0
    }
}

# Tidal periods (hours)
M2_PERIOD = 12.42  # Principal lunar semi-diurnal
K1_PERIOD = 23.93  # Lunar diurnal
S2_PERIOD = 12.00  # Principal solar semi-diurnal - NEW
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

def write_discharge_file(output_dir, filename, split_ratio):
    """
    Generate daily discharge time series for one year
    
    Args:
        output_dir: Output directory path
        filename: Output CSV filename
        split_ratio: Fraction of total discharge for this branch
    """
    filepath = os.path.join(output_dir, filename)
    print(f"Generating {filepath} (Flow split: {split_ratio*100:.0f}%)...")
    
    random.seed(RANDOM_SEED)
    
    with open(filepath, 'w', newline='') as f:
        # Write header (required by model)
        f.write('Date,Q\n')
        
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
            
            # Write to file - use actual newline
            f.write(f"{current_date.strftime('%Y-%m-%d')},{Q_branch:.2f}\n")
            
            current_date += timedelta(days=1)
    
    print(f"  ✓ Generated {filename} (365 days)")

# ============================================================================
# TIDAL GENERATION
# ============================================================================

def write_tidal_file(output_dir, filename, params):
    """
    Generate hourly tidal elevation time series for one year
    
    Uses harmonic analysis: H(t) = mean + M2*sin() + K1*sin() + S2*sin() + spring-neap
    
    PHASE 5 UPGRADE (Dec 2025): Added S2 constituent for proper spring-neap forcing.
    The South China Sea is a mixed semidiurnal system with form factor F ≈ 0.25-0.50.
    
    Reference: Nguyen et al. (2006) Continental Shelf Research
    
    Args:
        output_dir: Output directory path
        filename: Output CSV filename
        params: Dictionary with tidal constituent amplitudes and phases
    """
    filepath = os.path.join(output_dir, filename)
    print(f"Generating {filepath}...")
    
    random.seed(RANDOM_SEED)
    
    with open(filepath, 'w', newline='') as f:
        # Write header (required by model)
        f.write('DateTime,H\n')
        
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
            
            # S2 component (solar semi-diurnal, 12.00 hour period) - NEW
            # This creates proper spring-neap modulation via M2-S2 beat frequency
            A_S2 = params.get('A_S2', 0.0)
            phi_S2 = params.get('phi_S2', 0.0)
            if A_S2 > 0:
                h += A_S2 * math.sin(
                    2 * math.pi * hour_index / S2_PERIOD + phi_S2
                )
            
            # Spring-neap modulation (natural from M2+S2, but add residual)
            # Only apply if S2 is NOT already included (backward compatibility)
            if A_S2 == 0:
                spring_neap_factor = 1.0 + 0.1 * math.sin(
                    2 * math.pi * hour_index / SPRING_NEAP_PERIOD
                )
                h *= spring_neap_factor
            
            # Add small atmospheric/weather noise (±5 cm)
            h += random.gauss(0, 0.05)
            
            # Write to file - use actual newline
            f.write(f"{current_time.strftime('%Y-%m-%d %H:%M')},{h:.4f}\n")
            
            current_time += timedelta(hours=1)
            hour_index += 1
    
    print(f"  ✓ Generated {filename} (8760 hours)")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def generate_tien_river_forcing():
    """Generate all forcing files for Tien River case"""
    print("\n" + "-" * 70)
    print("Generating forcing for: Tien_River Case")
    print("-" * 70)
    
    ensure_directory(TIEN_FORCING_DIR)
    
    # Discharge
    write_discharge_file(TIEN_FORCING_DIR, 'Tien_Input.csv', RATIO_TIEN)
    
    # Tidal files
    for filename, params in TIEN_TIDAL_STATIONS.items():
        write_tidal_file(TIEN_FORCING_DIR, filename, params)

def generate_hau_river_forcing():
    """Generate all forcing files for Hau River case"""
    print("\n" + "-" * 70)
    print("Generating forcing for: Hau_River Case")
    print("-" * 70)
    
    ensure_directory(HAU_FORCING_DIR)
    
    # Discharge
    write_discharge_file(HAU_FORCING_DIR, 'Hau_Input.csv', RATIO_HAU)
    
    # Tidal files
    for filename, params in HAU_TIDAL_STATIONS.items():
        write_tidal_file(HAU_FORCING_DIR, filename, params)

def generate_root_forcing():
    """Generate forcing files in root forcing_data directory for convenience"""
    print("\n" + "-" * 70)
    print("Generating forcing in: root forcing_data (convenience copy)")
    print("-" * 70)
    
    ensure_directory(ROOT_FORCING_DIR)
    
    # Both discharge files
    write_discharge_file(ROOT_FORCING_DIR, 'Tien_Input.csv', RATIO_TIEN)
    write_discharge_file(ROOT_FORCING_DIR, 'Hau_Input.csv', RATIO_HAU)
    
    # All tidal files
    for filename, params in TIEN_TIDAL_STATIONS.items():
        write_tidal_file(ROOT_FORCING_DIR, filename, params)
    for filename, params in HAU_TIDAL_STATIONS.items():
        write_tidal_file(ROOT_FORCING_DIR, filename, params)

def main():
    """
    Main execution: Generate all forcing files for both cases
    """
    print("=" * 70)
    print("MEKONG DELTA FORCING GENERATOR")
    print("Parallel Estuary Approach - Savenije (2012)")
    print("=" * 70)
    
    # Generate forcing for each case
    generate_tien_river_forcing()
    generate_hau_river_forcing()
    generate_root_forcing()
    
    # Summary
    print("\n" + "=" * 70)
    print("GENERATION COMPLETE")
    print("=" * 70)
    print()
    print("Files generated for Tien_River case:")
    print(f"  Directory: {TIEN_FORCING_DIR}/")
    print("    - Tien_Input.csv (60% of total discharge)")
    print("    - MyTho_Tide.csv")
    print("    - HamLuong_Tide.csv")
    print("    - CoChien_Tide.csv")
    print()
    print("Files generated for Hau_River case:")
    print(f"  Directory: {HAU_FORCING_DIR}/")
    print("    - Hau_Input.csv (40% of total discharge)")
    print("    - Hau_Combined_Tide.csv")
    print()
    print("Discharge characteristics (2017):")
    print(f"  - Dry season base: {Q_BASE:.0f} m³/s")
    print(f"  - Flood peak: {Q_BASE + Q_FLOOD_PEAK:.0f} m³/s")
    print(f"  - Tien branch max: {(Q_BASE + Q_FLOOD_PEAK) * RATIO_TIEN:.0f} m³/s")
    print(f"  - Hau branch max: {(Q_BASE + Q_FLOOD_PEAK) * RATIO_HAU:.0f} m³/s")
    print()
    print("Tidal characteristics:")
    print(f"  - M2 period: {M2_PERIOD:.2f} hours")
    print(f"  - K1 period: {K1_PERIOD:.2f} hours")
    print("=" * 70)

if __name__ == '__main__':
    main()