#!/usr/bin/env python3
"""
Generate Lateral Nutrient/Carbon Loads from Land Use Map
=========================================================

This script converts the land use map (from satellite/GIS analysis) into
spatially-explicit lateral loads for the C-GEM biogeochemistry model.

EMISSION FACTORS (Literature-Based):
====================================

URBAN (Sewage + Stormwater):
- NH4: 150 kg N/ha/yr (Garnier et al., 2005 - Seine)
- NO3: 50 kg N/ha/yr 
- PO4: 20 kg P/ha/yr
- TOC: 500 kg C/ha/yr (BOD proxy)
- Source: IPCC 2006, Kroeze et al. (2002)

AQUACULTURE (Shrimp/Fish Ponds):
- NH4: 200 kg N/ha/yr (Excess feed + excretion)
- NO3: 100 kg N/ha/yr
- PO4: 40 kg P/ha/yr
- TOC: 1000 kg C/ha/yr (Uneaten feed, feces)
- Source: Páez-Osuna (2001), Naylor et al. (2000)

RICE PADDIES (Fertilizer Runoff):
- NH4: 30 kg N/ha/yr (Urea breakdown)
- NO3: 20 kg N/ha/yr
- PO4: 5 kg P/ha/yr
- TOC: 100 kg C/ha/yr (Straw decomposition)
- Source: Yan et al. (2003), MRC (2018)

MANGROVES (Tidal Carbon Export):
- TOC: 2000 kg C/ha/yr (Net tidal export of DOC)
- NH4: -10 kg N/ha/yr (N sink - denitrification)
- NO3: -20 kg N/ha/yr (N sink)
- PO4: 2 kg P/ha/yr
- Source: Alongi (2014), Bouillon et al. (2008)

FRUIT ORCHARDS:
- NH4: 40 kg N/ha/yr
- NO3: 30 kg N/ha/yr
- PO4: 10 kg P/ha/yr
- TOC: 200 kg C/ha/yr
- Source: Based on regional fertilizer application rates

UNIT CONVERSIONS:
=================
- kg/ha/yr → g/s per km² riparian zone
- 1 kg/ha/yr = 1000 g / (10000 m² × 365.25 × 24 × 3600 s)
            = 3.17e-9 g/m²/s
- For 1 km² = 1e6 m², multiply by 1e6 → 3.17e-3 g/s per km² per (kg/ha/yr)

LATERAL FLOW ESTIMATION:
========================
- Assume 0.001 m³/s/km² diffuse runoff (dry season)
- Wet season: 0.01 m³/s/km²
- This represents drainage canal + groundwater seepage

References:
- Garnier et al. (2005) Biogeochemistry 77, 213-242
- Alongi (2014) Annual Review Mar. Sci. 6, 195-219
- Páez-Osuna (2001) Environment International 27, 451-462
- MRC (2018) State of Basin Report

Author: CGEM Development Team
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

# Grid spacing (from topology)
DX_M = 2000.0  # Grid cell size [m]

# ===========================================================================
# EMISSION FACTORS [kg/ha/yr]
# ===========================================================================

# Format: {LandUse: {Species: emission_kg_ha_yr}}
EMISSION_FACTORS = {
    "Urban": {
        "NH4": 150.0,    # Sewage + stormwater
        "NO3": 50.0,
        "PO4": 20.0,
        "TOC": 500.0,    # BOD proxy
        "DIC": 200.0,    # Carbonate from concrete/roads
    },
    "Rice": {
        "NH4": 30.0,     # Urea breakdown
        "NO3": 20.0,     # Nitrification of fertilizer
        "PO4": 5.0,
        "TOC": 100.0,    # Straw/stubble decomposition
        "DIC": 50.0,
    },
    "Aqua": {
        "NH4": 200.0,    # Fish excretion + uneaten feed
        "NO3": 100.0,
        "PO4": 40.0,
        "TOC": 1000.0,   # High organic load
        "DIC": 300.0,
    },
    "Mangrove": {
        "NH4": -10.0,    # Net N sink (denitrification)
        "NO3": -20.0,    # Strong N sink
        "PO4": 2.0,      # Slight P release
        "TOC": 2000.0,   # Major carbon export (tidal)
        "DIC": 500.0,    # Respiration CO2
    },
    "Fruit": {
        "NH4": 40.0,
        "NO3": 30.0,
        "PO4": 10.0,
        "TOC": 200.0,
        "DIC": 100.0,
    },
}

# Conversion: kg/ha/yr to g/s per km²
# 1 kg/ha/yr × (1000 g/kg) × (1 ha/10000 m²) × (1 yr/31557600 s) × (1e6 m²/km²)
KG_HA_YR_TO_G_S_KM2 = 1000.0 / 10000.0 / 31557600.0 * 1e6  # ≈ 3.17e-3

# Lateral flow rate [m³/s per km²]
# Dry season: ~0.001 m³/s/km² (groundwater + drainage)
# Wet season: ~0.01 m³/s/km²
# Point sources (cities): ~0.05 m³/s/km² (sewage outfalls)
# 
# Using moderate rate (0.002) for realistic dry season
LATERAL_FLOW_RATE = 0.002  # m³/s per km² of riparian zone


def load_landuse_map() -> pd.DataFrame:
    """Load the land use map CSV."""
    path = CASE_DIR / "landuse_map.csv"
    if not path.exists():
        raise FileNotFoundError(f"Land use map not found: {path}\n"
                                "Run generate_synthetic_landuse.py first!")
    return pd.read_csv(path)


def calculate_loads(landuse_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate lateral loads from land use percentages.
    
    Returns DataFrame with columns:
    - Branch, Segment_Index, Distance_km
    - Q_lat_m3_s: Lateral inflow [m³/s]
    - NH4_load_g_s, NO3_load_g_s, PO4_load_g_s, TOC_load_g_s, DIC_load_g_s
    """
    records = []
    
    for _, row in landuse_df.iterrows():
        branch = row["Branch"]
        dist_km = row["Distance_km"]
        area_km2 = row["Segment_Area_km2"]
        
        # Calculate segment index (for C-code grid mapping)
        # Distance in km → grid index (assuming DX_M spacing)
        segment_idx = int(dist_km * 1000 / DX_M)
        
        # Lateral flow for this segment
        Q_lat = LATERAL_FLOW_RATE * area_km2
        
        # Calculate mass loads from each land use type
        loads = {sp: 0.0 for sp in ["NH4", "NO3", "PO4", "TOC", "DIC"]}
        
        for lu_type in ["Urban", "Rice", "Aqua", "Mangrove", "Fruit"]:
            pct = row.get(f"Pct_{lu_type}", 0.0) / 100.0  # Convert to fraction
            area_lu = area_km2 * pct  # km² of this land use
            
            for species, ef in EMISSION_FACTORS[lu_type].items():
                # Mass flux [g/s] = EF [kg/ha/yr] × conversion × area [km²]
                mass_flux = ef * KG_HA_YR_TO_G_S_KM2 * area_lu
                loads[species] += mass_flux
        
        # Calculate concentrations [g/m³ = mg/L]
        # This will be used for species that need concentration BC
        conc = {}
        if Q_lat > 1e-10:
            for sp, mass in loads.items():
                conc[sp] = mass / Q_lat  # g/s ÷ m³/s = g/m³
        else:
            for sp in loads:
                conc[sp] = 0.0
        
        records.append({
            "Branch": branch,
            "Segment_Index": segment_idx,
            "Distance_km": dist_km,
            "Area_km2": area_km2,
            "Q_lat_m3_s": Q_lat,
            "NH4_load_g_s": loads["NH4"],
            "NO3_load_g_s": loads["NO3"],
            "PO4_load_g_s": loads["PO4"],
            "TOC_load_g_s": loads["TOC"],
            "DIC_load_g_s": loads["DIC"],
            "NH4_conc_mg_L": conc["NH4"],
            "NO3_conc_mg_L": conc["NO3"],
            "PO4_conc_mg_L": conc["PO4"],
            "TOC_conc_mg_L": conc["TOC"],
            "DIC_conc_mg_L": conc["DIC"],
        })
    
    return pd.DataFrame(records)


def aggregate_to_grid(loads_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate 1-km segments to model grid cells (2 km spacing).
    
    The C-code needs loads at each grid cell, so we sum adjacent segments.
    """
    # Group by branch and segment index, sum the loads
    agg_df = loads_df.groupby(["Branch", "Segment_Index"]).agg({
        "Distance_km": "mean",
        "Area_km2": "sum",
        "Q_lat_m3_s": "sum",
        "NH4_load_g_s": "sum",
        "NO3_load_g_s": "sum",
        "PO4_load_g_s": "sum",
        "TOC_load_g_s": "sum",
        "DIC_load_g_s": "sum",
    }).reset_index()
    
    # Recalculate concentrations for aggregated cells
    for sp in ["NH4", "NO3", "PO4", "TOC", "DIC"]:
        agg_df[f"{sp}_conc_mg_L"] = np.where(
            agg_df["Q_lat_m3_s"] > 1e-10,
            agg_df[f"{sp}_load_g_s"] / agg_df["Q_lat_m3_s"],
            0.0
        )
    
    return agg_df


def main():
    print("=" * 70)
    print("GENERATING LATERAL LOADS FROM LAND USE MAP")
    print("=" * 70)
    
    # Load land use map
    print("\nLoading land use map...")
    landuse_df = load_landuse_map()
    print(f"  Loaded {len(landuse_df)} segments from {len(landuse_df['Branch'].unique())} branches")
    
    # Calculate loads
    print("\nCalculating lateral loads...")
    loads_df = calculate_loads(landuse_df)
    
    # Aggregate to grid
    print("Aggregating to model grid (2 km cells)...")
    grid_loads = aggregate_to_grid(loads_df)
    
    # Save output
    output_path = CASE_DIR / "lateral_sources.csv"
    grid_loads.to_csv(output_path, index=False)
    print(f"\nSaved lateral sources to: {output_path}")
    
    # Print summary
    print("\n" + "=" * 70)
    print("LATERAL LOAD SUMMARY BY BRANCH")
    print("=" * 70)
    
    summary = grid_loads.groupby("Branch").agg({
        "Q_lat_m3_s": "sum",
        "NH4_load_g_s": "sum",
        "NO3_load_g_s": "sum",
        "PO4_load_g_s": "sum",
        "TOC_load_g_s": "sum",
        "DIC_load_g_s": "sum",
    })
    
    # Convert to daily loads for readability
    summary_daily = summary.copy()
    summary_daily["Q_lat_m3_day"] = summary["Q_lat_m3_s"] * 86400
    for sp in ["NH4", "NO3", "PO4", "TOC", "DIC"]:
        summary_daily[f"{sp}_kg_day"] = summary[f"{sp}_load_g_s"] * 86400 / 1000
    
    print("\nDaily loads by branch:")
    print(summary_daily[["Q_lat_m3_day", "NH4_kg_day", "NO3_kg_day", 
                         "PO4_kg_day", "TOC_kg_day", "DIC_kg_day"]].round(1).to_string())
    
    # Identify hotspots
    print("\n" + "=" * 70)
    print("POLLUTION HOTSPOTS (NH4 load > 0.01 g/s)")
    print("=" * 70)
    
    hotspots = grid_loads[grid_loads["NH4_load_g_s"] > 0.01].sort_values(
        "NH4_load_g_s", ascending=False).head(10)
    
    if len(hotspots) > 0:
        for _, row in hotspots.iterrows():
            print(f"  {row['Branch']:15s} km {row['Distance_km']:5.1f}: "
                  f"NH4={row['NH4_load_g_s']*1000:.1f} mg/s, "
                  f"TOC={row['TOC_load_g_s']*1000:.1f} mg/s")
    
    # Mangrove carbon export
    print("\n" + "=" * 70)
    print("MANGROVE CARBON EXPORT (TOC load > 0.05 g/s)")
    print("=" * 70)
    
    mangrove_export = grid_loads[grid_loads["TOC_load_g_s"] > 0.05].sort_values(
        "TOC_load_g_s", ascending=False).head(10)
    
    if len(mangrove_export) > 0:
        for _, row in mangrove_export.iterrows():
            print(f"  {row['Branch']:15s} km {row['Distance_km']:5.1f}: "
                  f"TOC={row['TOC_load_g_s']*1000:.1f} mg/s")
    
    # Total delta loads
    print("\n" + "=" * 70)
    print("TOTAL MEKONG DELTA LATERAL LOADS")
    print("=" * 70)
    
    totals = summary.sum()
    print(f"  Total lateral flow:  {totals['Q_lat_m3_s']*86400:.0f} m³/day")
    print(f"  Total NH4 load:      {totals['NH4_load_g_s']*86400/1000:.1f} kg N/day")
    print(f"  Total NO3 load:      {totals['NO3_load_g_s']*86400/1000:.1f} kg N/day")
    print(f"  Total PO4 load:      {totals['PO4_load_g_s']*86400/1000:.1f} kg P/day")
    print(f"  Total TOC load:      {totals['TOC_load_g_s']*86400/1000:.1f} kg C/day")
    print(f"  Total DIC load:      {totals['DIC_load_g_s']*86400/1000:.1f} kg C/day")
    
    return grid_loads


if __name__ == "__main__":
    main()
