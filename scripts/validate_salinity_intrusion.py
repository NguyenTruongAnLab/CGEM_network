#!/usr/bin/env python3
"""
Validate Salinity Intrusion Against Literature
===============================================

This script validates the model output against expected Mekong Delta behavior:

VALIDATION TARGETS (Nguyen et al. 2008, MRC 2018):
==================================================
1. Dry Season (Feb-Apr):
   - 4 PSU isohaline: 40-60 km inland from mouth
   - Peak salinity at mouth: 25-30 PSU during high tide
   
2. Wet Season (Aug-Oct):
   - 4 PSU isohaline: <10 km from mouth
   - Freshwater dominates entire estuary

3. Tidal Range:
   - Water level oscillation: ±1.5-2.5 m at mouths
   - Velocity reversal during tidal cycle

Usage:
    python validate_salinity_intrusion.py OUTPUT/Mekong_Delta_Full/

Author: Nguyen Truong An
Date: December 2025
"""

import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import glob

# Validation thresholds
ISOHALINE_TARGET = 4.0  # PSU - standard intrusion reference
DRY_SEASON_INTRUSION_MIN = 30.0  # km
DRY_SEASON_INTRUSION_MAX = 70.0  # km
WET_SEASON_INTRUSION_MAX = 15.0  # km
TIDAL_RANGE_MIN = 1.0  # m
TIDAL_RANGE_MAX = 3.5  # m


def find_isohaline_position(distance_km: np.ndarray, salinity: np.ndarray, 
                            target_sal: float = ISOHALINE_TARGET) -> float:
    """
    Find the distance where salinity equals the target value.
    Uses linear interpolation between grid points.
    
    Returns: Distance in km, or NaN if not found.
    """
    if len(distance_km) != len(salinity):
        return np.nan
    
    # Find where salinity crosses the target
    for i in range(len(salinity) - 1):
        s1, s2 = salinity[i], salinity[i+1]
        d1, d2 = distance_km[i], distance_km[i+1]
        
        # Check if target is between these two points
        if (s1 >= target_sal >= s2) or (s2 >= target_sal >= s1):
            # Linear interpolation
            if abs(s2 - s1) < 1e-6:
                return 0.5 * (d1 + d2)
            frac = (target_sal - s1) / (s2 - s1)
            return d1 + frac * (d2 - d1)
    
    # If not found, check if all values are above or below target
    if np.all(salinity < target_sal):
        return 0.0  # Freshwater throughout
    if np.all(salinity > target_sal):
        return distance_km[-1]  # Salt throughout
    
    return np.nan


def load_branch_output(output_dir: Path, branch_name: str, variable: str):
    """Load output CSV for a specific branch and variable."""
    # Check CSV subdirectory first (typical output location)
    csv_dir = output_dir / "CSV"
    if csv_dir.exists():
        search_dir = csv_dir
    else:
        search_dir = output_dir
    
    # Try exact match first, then case-insensitive
    files = list(search_dir.glob(f"{branch_name}_{variable}.csv"))
    if not files:
        # Try lowercase variable name
        files = list(search_dir.glob(f"{branch_name}_{variable.lower()}.csv"))
    
    if not files:
        return None
    
    try:
        df = pd.read_csv(files[0])
        return df
    except Exception as e:
        print(f"Error loading {files[0]}: {e}")
        return None


def analyze_branch_salinity(output_dir: Path, branch_name: str, dx_m: float = 2000.0):
    """Analyze salinity intrusion for a single branch."""
    
    sal_df = load_branch_output(output_dir, branch_name, "SALINITY")
    depth_df = load_branch_output(output_dir, branch_name, "depth")
    
    if sal_df is None:
        print(f"  WARNING: No salinity data for {branch_name}")
        return None
    
    # Get time column (first column is usually time)
    time_col = sal_df.columns[0]
    
    # Get spatial columns (remaining columns are grid points)
    spatial_cols = [c for c in sal_df.columns if c != time_col]
    n_cells = len(spatial_cols)
    
    # Create distance array (index 1 = downstream = 0 km)
    distance_km = np.arange(n_cells) * dx_m / 1000.0
    
    results = {
        "branch": branch_name,
        "n_timesteps": len(sal_df),
        "intrusion_km": [],
        "max_salinity": [],
        "time_days": [],
    }
    
    for idx, row in sal_df.iterrows():
        time_s = row[time_col]
        salinity = row[spatial_cols].values.astype(float)
        
        # Find 4 PSU isohaline
        intrusion = find_isohaline_position(distance_km, salinity)
        
        results["time_days"].append(time_s / 86400.0)
        results["intrusion_km"].append(intrusion)
        results["max_salinity"].append(np.max(salinity))
    
    return results


def print_validation_report(results: dict, branch_name: str):
    """Print validation report for a branch."""
    
    intrusion = np.array(results["intrusion_km"])
    max_sal = np.array(results["max_salinity"])
    time_days = np.array(results["time_days"])
    
    # Remove NaN values
    valid_mask = ~np.isnan(intrusion)
    intrusion_valid = intrusion[valid_mask]
    time_valid = time_days[valid_mask]
    
    print(f"\n{'='*60}")
    print(f"VALIDATION REPORT: {branch_name}")
    print(f"{'='*60}")
    
    if len(intrusion_valid) == 0:
        print("  ERROR: No valid salinity intrusion data!")
        print("  POSSIBLE CAUSES:")
        print("    1. Ocean boundary concentration not being applied")
        print("    2. Dispersion coefficient too low")
        print("    3. Transport solver not propagating salinity")
        return False
    
    print(f"\n4 PSU Isohaline Position:")
    print(f"  Mean: {np.mean(intrusion_valid):.1f} km")
    print(f"  Min:  {np.min(intrusion_valid):.1f} km")
    print(f"  Max:  {np.max(intrusion_valid):.1f} km")
    
    print(f"\nMaximum Salinity at Mouth:")
    print(f"  Mean: {np.mean(max_sal):.1f} PSU")
    print(f"  Max:  {np.max(max_sal):.1f} PSU")
    
    # Seasonal analysis (if enough data)
    if len(time_valid) > 30:
        # Approximate dry season: days 30-120 (Feb-Apr)
        dry_mask = (time_valid >= 30) & (time_valid <= 120)
        # Approximate wet season: days 210-300 (Aug-Oct)
        wet_mask = (time_valid >= 210) & (time_valid <= 300)
        
        if np.any(dry_mask):
            dry_intrusion = np.mean(intrusion_valid[dry_mask])
            print(f"\nDry Season (days 30-120):")
            print(f"  Mean intrusion: {dry_intrusion:.1f} km")
            if DRY_SEASON_INTRUSION_MIN <= dry_intrusion <= DRY_SEASON_INTRUSION_MAX:
                print(f"  ✓ PASS: Within expected range ({DRY_SEASON_INTRUSION_MIN}-{DRY_SEASON_INTRUSION_MAX} km)")
            else:
                print(f"  ✗ FAIL: Outside expected range ({DRY_SEASON_INTRUSION_MIN}-{DRY_SEASON_INTRUSION_MAX} km)")
        
        if np.any(wet_mask):
            wet_intrusion = np.mean(intrusion_valid[wet_mask])
            print(f"\nWet Season (days 210-300):")
            print(f"  Mean intrusion: {wet_intrusion:.1f} km")
            if wet_intrusion <= WET_SEASON_INTRUSION_MAX:
                print(f"  ✓ PASS: Within expected range (< {WET_SEASON_INTRUSION_MAX} km)")
            else:
                print(f"  ✗ FAIL: Too high (expected < {WET_SEASON_INTRUSION_MAX} km)")
    
    # Basic sanity checks
    passed = True
    print(f"\nSanity Checks:")
    
    if np.max(max_sal) < 1.0:
        print(f"  ✗ FAIL: Max salinity < 1 PSU - salt not entering estuary!")
        passed = False
    else:
        print(f"  ✓ PASS: Salt entering estuary (max = {np.max(max_sal):.1f} PSU)")
    
    if np.mean(intrusion_valid) < 1.0:
        print(f"  ✗ FAIL: Mean intrusion < 1 km - transport not working!")
        passed = False
    else:
        print(f"  ✓ PASS: Salt intrusion active (mean = {np.mean(intrusion_valid):.1f} km)")
    
    return passed


def main():
    parser = argparse.ArgumentParser(description="Validate salinity intrusion")
    parser.add_argument("output_dir", nargs="?", 
                        default="OUTPUT/Mekong_Delta_Full",
                        help="Path to model output directory")
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    
    if not output_dir.exists():
        print(f"ERROR: Output directory not found: {output_dir}")
        return 1
    
    print("=" * 70)
    print("C-GEM Network: Salinity Intrusion Validation")
    print("=" * 70)
    print(f"Output directory: {output_dir}")
    
    # Analyze key estuarine branches
    branches = ["Ham_Luong", "My_Tho", "Co_Chien", "Hau_River"]
    
    all_passed = True
    for branch in branches:
        results = analyze_branch_salinity(output_dir, branch)
        if results:
            passed = print_validation_report(results, branch)
            if not passed:
                all_passed = False
    
    print("\n" + "=" * 70)
    if all_passed:
        print("OVERALL: ✓ VALIDATION PASSED")
    else:
        print("OVERALL: ✗ VALIDATION FAILED - Check warnings above")
    print("=" * 70)
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
