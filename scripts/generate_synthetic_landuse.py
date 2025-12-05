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
Date: November 2025
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

# Assumed riparian corridor width (km) for load calculation
RIPARIAN_WIDTH_KM = 2.0  # 1 km each side


def read_branch_lengths_from_topology(case_dir: Path) -> dict:
    """
    Read branch lengths from topology.csv file.
    
    Returns dict mapping branch name to length in km.
    """
    topology_path = case_dir / "topology.csv"
    if not topology_path.exists():
        raise FileNotFoundError(f"Topology file not found: {topology_path}")
    
    branch_lengths = {}
    with open(topology_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue
            
            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 5:
                try:
                    # Format: ID, Name, NodeUp, NodeDown, Length_m, ...
                    branch_name = parts[1]
                    length_m = float(parts[4])
                    length_km = int(length_m / 1000)
                    branch_lengths[branch_name] = length_km
                except (ValueError, IndexError):
                    continue
    
    return branch_lengths

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
    # Total length: 166km (from ocean mouth to Vam Nao junction)
    # 
    # DISTANCE CONVENTION: km 0 = OCEAN MOUTH, km 166 = UPSTREAM (Vam Nao)
    # 
    # GEOGRAPHIC ZONES (Literature: Nguyen et al. 2008, MRC 2018):
    # km 0-30: COASTAL - Dinh An/Tran De mouth, mangroves, shrimp farms
    # km 30-60: LOWER HAU - Soc Trang province, aquaculture/rice transition
    # km 60-100: CAN THO AREA - metropolitan zone (city center ~km 80 from ocean)
    #            Can Tho: Vietnam's 4th largest city, 1.5M population
    # km 100-140: MID HAU - Peri-urban, mixed agriculture
    # km 140-166: UPPER HAU - Near Vam Nao, rural agriculture (Long Xuyen area)
    if branch == "Hau_River":
        if distance_km < 30:
            # COASTAL - Dinh An/Tran De mouth, mangroves and shrimp farms
            lu = {"Pct_Urban": 5, "Pct_Rice": 10, "Pct_Aqua": 45, 
                  "Pct_Mangrove": 35, "Pct_Fruit": 5}
        elif 30 <= distance_km < 60:
            # LOWER HAU - Soc Trang, aquaculture transition
            lu = {"Pct_Urban": 15, "Pct_Rice": 25, "Pct_Aqua": 40, 
                  "Pct_Mangrove": 15, "Pct_Fruit": 5}
        elif 60 <= distance_km < 75:
            # Approaching Can Tho - peri-urban/industrial
            lu = {"Pct_Urban": 40, "Pct_Rice": 30, "Pct_Aqua": 15, 
                  "Pct_Mangrove": 5, "Pct_Fruit": 10}
        elif 75 <= distance_km <= 95:
            # CAN THO metropolitan core (city center at ~km 85 from ocean)
            lu = {"Pct_Urban": 70, "Pct_Rice": 15, "Pct_Aqua": 5, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 10}
        elif 95 < distance_km <= 120:
            # Upstream of Can Tho - peri-urban agriculture
            lu = {"Pct_Urban": 35, "Pct_Rice": 40, "Pct_Aqua": 5, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 20}
        elif 120 < distance_km <= 145:
            # MID-UPPER HAU - rural agriculture
            lu = {"Pct_Urban": 15, "Pct_Rice": 70, "Pct_Aqua": 0, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 15}
        else:
            # UPPER HAU (km > 145) - near Vam Nao, rural/Long Xuyen
            lu = {"Pct_Urban": 10, "Pct_Rice": 75, "Pct_Aqua": 0, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 15}
    
    # ==== MY THO (Northernmost Tien - My Tho City) ====
    # Total length: 72km (from ocean mouth to Tien_Connector junction)
    #
    # DISTANCE CONVENTION: km 0 = OCEAN MOUTH, km 72 = UPSTREAM (junction)
    #
    # GEOGRAPHIC ZONES:
    # km 0-20: COASTAL - Cua Tieu/Cua Dai mouth, mangroves and aquaculture
    # km 20-40: My Tho city (Tien Giang capital, ~500k population, centered ~km 30)
    # km 40-55: Mid-reach fruit orchards
    # km 55-72: UPSTREAM - near junction, fruit orchards and rice
    elif branch == "My_Tho":
        if distance_km < 20:
            # COASTAL - mangroves and aquaculture
            lu = {"Pct_Urban": 10, "Pct_Rice": 15, "Pct_Aqua": 35, 
                  "Pct_Mangrove": 30, "Pct_Fruit": 10}
        elif 20 <= distance_km < 40:
            # My Tho city area (~500k population, city center ~km 30)
            lu = {"Pct_Urban": 55, "Pct_Rice": 20, "Pct_Aqua": 5, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 20}
        elif 40 <= distance_km < 55:
            # Mid-reach fruit orchards
            lu = {"Pct_Urban": 10, "Pct_Rice": 30, "Pct_Aqua": 0, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 60}
        else:
            # UPSTREAM (km > 55) - near junction
            lu = {"Pct_Urban": 15, "Pct_Rice": 50, "Pct_Aqua": 0, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 35}
    
    # ==== HAM LUONG (Central - Ben Tre fruit region) ====
    # Total length: 76km (from ocean mouth to Tien_Connector junction)
    #
    # DISTANCE CONVENTION: km 0 = OCEAN MOUTH, km 76 = UPSTREAM (junction)
    #
    # GEOGRAPHIC ZONES:
    # km 0-25: COASTAL - Ham Luong mouth, mangroves and aquaculture  
    # km 25-45: Ben Tre city area (~200k population, city center ~km 35)
    # km 45-76: UPSTREAM - famous fruit orchard region (Ben Tre - "coconut land")
    elif branch == "Ham_Luong":
        if distance_km < 25:
            # COASTAL - mangroves and aquaculture
            lu = {"Pct_Urban": 5, "Pct_Rice": 10, "Pct_Aqua": 40, 
                  "Pct_Mangrove": 40, "Pct_Fruit": 5}
        elif 25 <= distance_km < 45:
            # Ben Tre city area (~200k population, city center ~km 35)
            lu = {"Pct_Urban": 40, "Pct_Rice": 20, "Pct_Aqua": 10, 
                  "Pct_Mangrove": 5, "Pct_Fruit": 25}
        else:
            # UPSTREAM (km > 45) - fruit orchard region
            lu = {"Pct_Urban": 5, "Pct_Rice": 20, "Pct_Aqua": 0, 
                  "Pct_Mangrove": 0, "Pct_Fruit": 75}
    
    # ==== CO CHIEN (South Tien - aquaculture zone) ====
    # Total length: 98km (from ocean mouth to Co_Chien_Split junction)
    #
    # DISTANCE CONVENTION: km 0 = OCEAN MOUTH, km 98 = UPSTREAM (junction)
    #
    # GEOGRAPHIC ZONES:
    # km 0-25: COASTAL - Cung Hau/Co Chien mouth, mangroves and aquaculture
    # km 25-50: Transition zone - aquaculture/agriculture mix
    # km 50-75: Vinh Long area (city center ~km 60, ~150k population)
    # km 75-98: UPSTREAM - near junction, mixed agriculture
    elif branch == "Co_Chien":
        if distance_km < 25:
            # COASTAL - heavy aquaculture and mangroves
            lu = {"Pct_Urban": 5, "Pct_Rice": 5, "Pct_Aqua": 55, 
                  "Pct_Mangrove": 30, "Pct_Fruit": 5}
        elif 25 <= distance_km < 50:
            # Transition zone
            lu = {"Pct_Urban": 10, "Pct_Rice": 25, "Pct_Aqua": 35, 
                  "Pct_Mangrove": 20, "Pct_Fruit": 10}
        elif 50 <= distance_km < 75:
            # Vinh Long area (city center ~km 60)
            lu = {"Pct_Urban": 25, "Pct_Rice": 35, "Pct_Aqua": 15, 
                  "Pct_Mangrove": 5, "Pct_Fruit": 20}
        else:
            # UPSTREAM (km > 75) - near junction, agriculture
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


def generate_landuse_map(case_dir: Path = CASE_DIR) -> pd.DataFrame:
    """
    Generate complete land use map for all branches.
    Reads branch lengths from topology.csv.
    """
    # Read branch lengths from topology file
    branch_lengths = read_branch_lengths_from_topology(case_dir)
    print(f"  Read {len(branch_lengths)} branches from topology.csv")
    for branch, length in branch_lengths.items():
        print(f"    {branch}: {length} km")
    
    records = []
    
    for branch, length_km in branch_lengths.items():
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
