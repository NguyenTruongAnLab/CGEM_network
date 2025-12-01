#!/usr/bin/env python3
"""
Generate Synthetic Land Use Map for Mekong Delta
=================================================

This script mimics the output of a GIS analysis on JAXA ALOS-2/PALSAR or
Sentinel-2 satellite data. It creates a CSV file mapping river kilometers 
to land use types for each branch.

LAND USE CLASSIFICATION (Simplified from JAXA Land Use/Cover Map):
==================================================================
1. URBAN       - Built-up areas, cities (Can Tho, My Tho, Ben Tre)
2. RICE        - Paddy fields (dominant in inland delta)
3. AQUACULTURE - Shrimp/fish ponds (coastal areas)
4. MANGROVE    - Coastal mangrove forests (carbon sink/source)
5. FRUIT       - Orchards (mango, longan, durian - middle delta)

SPATIAL PATTERNS (Based on Mekong Delta Geography):
===================================================
- Upstream (km 0-30):     Inland rice paddies, some fruit orchards
- Mid-reach (km 30-60):   Urban centers (Can Tho, My Tho), mixed agriculture
- Downstream (km 60+):    Coastal - mangroves, aquaculture

BRANCH-SPECIFIC FEATURES:
- Hau_River: Can Tho city at km 30-50 (major pollution source)
- My_Tho: My Tho city at km 40-55
- Ham_Luong: Ben Tre area, fruit orchards
- Co_Chien: Mixed rice/aquaculture
- Tien_Main/Lower: Inland agriculture

References:
- JAXA Earth Observation (https://www.eorc.jaxa.jp/ALOS/en/dataset/lulc_e.htm)
- Tuan et al. (2020) - Mekong Delta land use classification
- MRC State of Basin Report (2018)

Author: Nguyen Truong An
Date: November 2024
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np

# ===========================================================================
# CONFIGURATION
# ===========================================================================

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent
CASE_DIR = PROJECT_ROOT / "INPUT" / "Cases" / "Mekong_Delta_Full"

# Branch lengths (km) - from topology
BRANCH_LENGTHS = {
    "Tien_Main": 50,
    "Hau_Main": 45,
    "Vam_Nao": 7,
    "Tien_Lower": 50,
    "Tien_Connector": 15,
    "Co_Chien": 85,
    "My_Tho": 55,
    "Ham_Luong": 60,
    "Hau_River": 90,
}

# Assumed riparian corridor width (km) for load calculation
RIPARIAN_WIDTH_KM = 2.0  # 1 km each side

# ===========================================================================
# LAND USE ASSIGNMENT RULES
# ===========================================================================

def assign_landuse(branch: str, distance_km: float) -> dict:
    """
    Assign land use percentages based on branch and distance from upstream.
    
    Returns dict with: Pct_Urban, Pct_Rice, Pct_Aqua, Pct_Mangrove, Pct_Fruit
    All percentages should sum to 100.
    """
    
    # Default: Inland rice-dominated
    lu = {
        "Pct_Urban": 5,
        "Pct_Rice": 85,
        "Pct_Aqua": 0,
        "Pct_Mangrove": 0,
        "Pct_Fruit": 10,
    }
    
    # ==== HAU RIVER (Main Hau - Can Tho City) ====
    if branch == "Hau_River":
        if 25 <= distance_km <= 50:
            # Can Tho metropolitan area (major city, 1.2M population)
            lu = {"Pct_Urban": 70, "Pct_Rice": 15, "Pct_Aqua": 5, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 10}
        elif 50 < distance_km <= 70:
            # Peri-urban / mixed
            lu = {"Pct_Urban": 30, "Pct_Rice": 40, "Pct_Aqua": 15, 
                  "Pct_Mangrove": 5, "Pct_Fruit": 10}
        elif distance_km > 70:
            # Coastal zone - aquaculture and mangroves
            lu = {"Pct_Urban": 5, "Pct_Rice": 10, "Pct_Aqua": 45, 
                  "Pct_Mangrove": 35, "Pct_Fruit": 5}
        else:
            # Upper Hau - agriculture
            lu = {"Pct_Urban": 10, "Pct_Rice": 75, "Pct_Aqua": 0, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 15}
    
    # ==== MY THO (Northernmost Tien - My Tho City) ====
    elif branch == "My_Tho":
        if 35 <= distance_km <= 50:
            # My Tho city area (~500k population)
            lu = {"Pct_Urban": 55, "Pct_Rice": 20, "Pct_Aqua": 5, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 20}
        elif distance_km > 50:
            # Coastal
            lu = {"Pct_Urban": 10, "Pct_Rice": 15, "Pct_Aqua": 35, 
                  "Pct_Mangrove": 30, "Pct_Fruit": 10}
        elif distance_km < 15:
            # Upper reach - near junction
            lu = {"Pct_Urban": 15, "Pct_Rice": 60, "Pct_Aqua": 0, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 25}
        else:
            # Mid-reach fruit orchards
            lu = {"Pct_Urban": 10, "Pct_Rice": 40, "Pct_Aqua": 0, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 50}
    
    # ==== HAM LUONG (Central - Ben Tre fruit region) ====
    elif branch == "Ham_Luong":
        if 40 <= distance_km <= 55:
            # Ben Tre city area (~200k population)
            lu = {"Pct_Urban": 35, "Pct_Rice": 25, "Pct_Aqua": 10, 
                  "Pct_Mangrove": 5, "Pct_Fruit": 25}
        elif distance_km > 55:
            # Coastal
            lu = {"Pct_Urban": 5, "Pct_Rice": 10, "Pct_Aqua": 40, 
                  "Pct_Mangrove": 40, "Pct_Fruit": 5}
        else:
            # Famous fruit orchard region (Ben Tre province)
            lu = {"Pct_Urban": 5, "Pct_Rice": 25, "Pct_Aqua": 0, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 70}
    
    # ==== CO CHIEN (South Tien - aquaculture zone) ====
    elif branch == "Co_Chien":
        if distance_km > 70:
            # Coastal - heavy aquaculture
            lu = {"Pct_Urban": 5, "Pct_Rice": 5, "Pct_Aqua": 55, 
                  "Pct_Mangrove": 30, "Pct_Fruit": 5}
        elif 50 < distance_km <= 70:
            # Transition zone
            lu = {"Pct_Urban": 10, "Pct_Rice": 30, "Pct_Aqua": 30, 
                  "Pct_Mangrove": 10, "Pct_Fruit": 20}
        else:
            # Inland rice/fruit
            lu = {"Pct_Urban": 10, "Pct_Rice": 55, "Pct_Aqua": 5, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 30}
    
    # ==== TIEN MAIN (Upper Tien - agriculture) ====
    elif branch == "Tien_Main":
        # Primarily agricultural (upstream of delta)
        lu = {"Pct_Urban": 10, "Pct_Rice": 80, "Pct_Aqua": 0, 
              "Pct_Mangrove": 0, "Pct_Fruit": 10}
    
    # ==== HAU MAIN (Upper Hau) ====
    elif branch == "Hau_Main":
        if distance_km > 30:
            # Near Long Xuyen city
            lu = {"Pct_Urban": 25, "Pct_Rice": 65, "Pct_Aqua": 0, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 10}
        else:
            lu = {"Pct_Urban": 5, "Pct_Rice": 85, "Pct_Aqua": 0, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 10}
    
    # ==== TIEN LOWER (Between Vam Nao and Co Chien split) ====
    elif branch == "Tien_Lower":
        lu = {"Pct_Urban": 10, "Pct_Rice": 70, "Pct_Aqua": 0, 
              "Pct_Mangrove": 0, "Pct_Fruit": 20}
    
    # ==== TIEN CONNECTOR (Short reach) ====
    elif branch == "Tien_Connector":
        lu = {"Pct_Urban": 15, "Pct_Rice": 60, "Pct_Aqua": 0, 
              "Pct_Mangrove": 0, "Pct_Fruit": 25}
    
    # ==== VAM NAO (Hydraulic connector - minimal land use) ====
    elif branch == "Vam_Nao":
        # Very short channel, mostly water
        lu = {"Pct_Urban": 5, "Pct_Rice": 90, "Pct_Aqua": 0, 
              "Pct_Mangrove": 0, "Pct_Fruit": 5}
    
    return lu


def generate_landuse_map() -> pd.DataFrame:
    """
    Generate complete land use map for all branches.
    """
    records = []
    
    for branch, length_km in BRANCH_LENGTHS.items():
        # Generate 1-km segments
        for dist_km in range(0, length_km + 1, 1):
            lu = assign_landuse(branch, dist_km)
            
            # Calculate area for this segment (km²)
            segment_area_km2 = 1.0 * RIPARIAN_WIDTH_KM  # 1 km segment × corridor width
            
            records.append({
                "Branch": branch,
                "Distance_km": dist_km,
                "Segment_Area_km2": segment_area_km2,
                **lu
            })
    
    return pd.DataFrame(records)


def main():
    print("=" * 70)
    print("GENERATING SYNTHETIC LAND USE MAP FOR MEKONG DELTA")
    print("=" * 70)
    
    # Create output directory
    CASE_DIR.mkdir(parents=True, exist_ok=True)
    
    # Generate land use map
    df = generate_landuse_map()
    
    # Save to CSV
    output_path = CASE_DIR / "landuse_map.csv"
    df.to_csv(output_path, index=False)
    print(f"\nSaved land use map to: {output_path}")
    
    # Print summary statistics
    print("\n" + "=" * 70)
    print("LAND USE SUMMARY BY BRANCH")
    print("=" * 70)
    
    summary = df.groupby("Branch").agg({
        "Pct_Urban": "mean",
        "Pct_Rice": "mean",
        "Pct_Aqua": "mean",
        "Pct_Mangrove": "mean",
        "Pct_Fruit": "mean",
        "Segment_Area_km2": "sum"
    }).round(1)
    
    print(summary.to_string())
    
    # Identify hotspots
    print("\n" + "=" * 70)
    print("POLLUTION HOTSPOTS (Urban > 50%)")
    print("=" * 70)
    
    hotspots = df[df["Pct_Urban"] > 50]
    if len(hotspots) > 0:
        for _, row in hotspots.iterrows():
            print(f"  {row['Branch']} km {row['Distance_km']}: {row['Pct_Urban']:.0f}% Urban")
    else:
        print("  None identified")
    
    print("\n" + "=" * 70)
    print("COASTAL ZONES (Mangrove > 20%)")
    print("=" * 70)
    
    coastal = df[df["Pct_Mangrove"] > 20]
    if len(coastal) > 0:
        branches_coastal = coastal.groupby("Branch")["Distance_km"].agg(["min", "max"])
        for branch, row in branches_coastal.iterrows():
            print(f"  {branch}: km {row['min']:.0f} - {row['max']:.0f}")
    
    return df


if __name__ == "__main__":
    main()
