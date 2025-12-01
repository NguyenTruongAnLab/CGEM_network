#!/usr/bin/env python3
"""
Generate Lateral Nutrient/Carbon Loads from Land Use Map
=========================================================

This script converts the land use map (from satellite/GIS analysis) into
spatially-explicit lateral loads for the C-GEM biogeochemistry model.
=======================================

1. VIRTUAL POLDERS (80/20 Solution):
   Instead of implementing complex gate logic in C code, this script generates
   time-varying lateral loads that mimic sluice gate operation:
   - Gate CLOSED (high tide): Q_lat = 0 (no pollution release)
   - Gate OPEN (low tide): Q_lat = HIGH (flushing pulse)
   
   This creates realistic "pollution pulses" without touching hydrodynamics.
   
2. SEASONAL HYDROLOGY:
   Runoff varies dramatically between seasons in the Mekong:
   - Dry season (Dec-May): 0.001 m³/s/km² (groundwater seepage only)
   - Wet season (Jun-Nov): 0.01 m³/s/km² (monsoon runoff)
   - Transition: Gradual interpolation

3. POINT SOURCES:
   Major cities (Can Tho, My Tho) treated as concentrated point sources
   with higher flow rates than diffuse runoff.

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
- Savenije (2012) Salinity and Tides in Alluvial Estuaries

Author: Nguyen Truong An
Date: December 2025 (Updated with Virtual Polders and Seasonal Hydrology)
"""

import argparse
import os
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Optional, Tuple

# ===========================================================================
# CONFIGURATION
# ===========================================================================

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent
CASE_DIR = PROJECT_ROOT / "INPUT" / "Cases" / "Mekong_Delta_Full"

# Grid spacing (from topology)
DX_M = 2000.0  # Grid cell size [m]

# Full year duration for comprehensive seasonal coverage
FULL_YEAR_DAYS = 365.0

# ===========================================================================
# SEASONAL HYDROLOGY PARAMETERS
# ===========================================================================

# Lateral flow rates [m³/s per km²] by season
SEASONAL_FLOW_RATES = {
    "dry": 0.001,       # Dec-May: minimal groundwater seepage
    "wet": 0.010,       # Jun-Nov: monsoon runoff
    "transition": 0.004, # Shoulder seasons
}

# Monthly flow rate multipliers (relative to dry season base)
# Based on Mekong River Commission data for delta rainfall/runoff
MONTHLY_FLOW_MULTIPLIERS = {
    1: 1.0,    # January - Dry
    2: 1.0,    # February - Dry
    3: 1.2,    # March - Late dry (early transition)
    4: 2.0,    # April - Transition
    5: 4.0,    # May - Early wet
    6: 8.0,    # June - Wet (monsoon onset)
    7: 10.0,   # July - Peak wet
    8: 10.0,   # August - Peak wet
    9: 9.0,    # September - Wet (flood peak)
    10: 6.0,   # October - Late wet
    11: 3.0,   # November - Transition
    12: 1.5,   # December - Early dry
}

# Dry/Wet season months (Mekong climate)
DRY_MONTHS = [12, 1, 2, 3, 4, 5]   # December to May
WET_MONTHS = [6, 7, 8, 9, 10, 11]  # June to November

# ===========================================================================
# VIRTUAL POLDER PARAMETERS (Gate Simulation)
# ===========================================================================

# Tidal threshold for gate operation [m above MSL]
# Gates typically close when water level exceeds ~0.5-1.0 m above mean
POLDER_TIDE_THRESHOLD = 0.5  # meters

# Flow multiplier when gate is OPEN (flushing)
# During low tide, accumulated waste is flushed out at higher rate
POLDER_FLUSH_MULTIPLIER = 3.0

# Time step for time-varying output [seconds]
# M2 tidal period = 12.42 hours = 44712 seconds
# Use 1-hour resolution for gate switching
OUTPUT_TIME_STEP = 3600.0  # 1 hour

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

# Point source cities (concentrated outfalls)
POINT_SOURCES = {
    "Can_Tho": {
        "branch": "Hau_Main",
        "distance_km": 80.0,  # Approximate distance from mouth
        "population": 1_500_000,
        "per_capita_NH4_g_day": 12.0,  # ~12g N/person/day
        "per_capita_TOC_g_day": 50.0,  # ~50g C/person/day (BOD)
    },
    "My_Tho": {
        "branch": "Tien_Main",
        "distance_km": 70.0,
        "population": 500_000,
        "per_capita_NH4_g_day": 12.0,
        "per_capita_TOC_g_day": 50.0,
    },
}

# Conversion: kg/ha/yr to g/s per km²
# 1 kg/ha/yr × (1000 g/kg) × (1 ha/10000 m²) × (1 yr/31557600 s) × (1e6 m²/km²)
KG_HA_YR_TO_G_S_KM2 = 1000.0 / 10000.0 / 31557600.0 * 1e6  # ≈ 3.17e-3


def load_landuse_map() -> pd.DataFrame:
    """Load the land use map CSV."""
    path = CASE_DIR / "landuse_map.csv"
    if not path.exists():
        raise FileNotFoundError(f"Land use map not found: {path}\n"
                                "Run generate_synthetic_landuse.py first!")
    return pd.read_csv(path)


def get_seasonal_flow_rate(season: str) -> float:
    """
    Get lateral flow rate based on season.
    
    Args:
        season: 'dry', 'wet', or 'transition'
    
    Returns:
        Flow rate [m³/s per km²]
    """
    if season not in SEASONAL_FLOW_RATES:
        print(f"Warning: Unknown season '{season}', using dry season rate")
        return SEASONAL_FLOW_RATES["dry"]
    return SEASONAL_FLOW_RATES[season]


def generate_tidal_signal(duration_days: float, dt: float = OUTPUT_TIME_STEP) -> np.ndarray:
    """
    Generate a synthetic M2 tidal water level signal.
    
    Args:
        duration_days: Simulation duration [days]
        dt: Time step [seconds]
    
    Returns:
        Array of water levels [m] at each timestep
    """
    # M2 tidal period = 12.42 hours
    M2_PERIOD_S = 12.42 * 3600.0
    omega = 2.0 * np.pi / M2_PERIOD_S
    
    # Generate time array
    n_steps = int(duration_days * 86400.0 / dt)
    time_s = np.arange(n_steps) * dt
    
    # M2 tide with ~1.5m amplitude (typical Mekong)
    # Add spring-neap modulation (14.76 day cycle)
    M2_amplitude = 1.5  # meters
    spring_neap_period = 14.76 * 86400.0  # seconds
    spring_neap_mod = 0.3 * np.sin(2 * np.pi * time_s / spring_neap_period)
    
    water_level = (M2_amplitude + spring_neap_mod) * np.sin(omega * time_s)
    
    return water_level


def calculate_polder_gate_factor(water_level: np.ndarray) -> np.ndarray:
    """
    Calculate gate operation factor based on water level.
    
    Virtual Polder Logic:
    - High tide (level > threshold): Gate CLOSED → factor = 0
    - Low tide (level < threshold): Gate OPEN → factor = FLUSH_MULTIPLIER
    - Smooth transition to avoid numerical discontinuities
    
    Args:
        water_level: Array of water levels [m]
    
    Returns:
        Array of flow multipliers (0 to FLUSH_MULTIPLIER)
    """
    # Smooth sigmoid transition (avoids step function)
    # Transition width ~0.2m around threshold
    transition_width = 0.2
    
    # Sigmoid: 0 when level > threshold, 1 when level < threshold
    gate_open = 1.0 / (1.0 + np.exp((water_level - POLDER_TIDE_THRESHOLD) / transition_width))
    
    # Scale by flush multiplier
    return gate_open * POLDER_FLUSH_MULTIPLIER


def calculate_loads(landuse_df: pd.DataFrame, season: str = "dry",
                    enable_polders: bool = False) -> pd.DataFrame:
    """
    Calculate lateral loads from land use percentages.
    
    Args:
        landuse_df: Land use DataFrame with Pct_* columns
        season: 'dry', 'wet', or 'transition'
        enable_polders: If True, output will support time-varying loads
    
    Returns DataFrame with columns:
    - Branch, Segment_Index, Distance_km
    - Q_lat_m3_s: Lateral inflow [m³/s]
    - NH4_load_g_s, NO3_load_g_s, PO4_load_g_s, TOC_load_g_s, DIC_load_g_s
    """
    records = []
    
    # Get season-appropriate flow rate
    base_flow_rate = get_seasonal_flow_rate(season)
    
    for _, row in landuse_df.iterrows():
        branch = row["Branch"]
        dist_km = row["Distance_km"]
        area_km2 = row["Segment_Area_km2"]
        
        # Calculate segment index (for C-code grid mapping)
        # Distance in km → grid index (assuming DX_M spacing)
        segment_idx = int(dist_km * 1000 / DX_M)
        
        # Check if this segment is in a "polder zone" (rice/aqua dominated)
        pct_polder = row.get("Pct_Rice", 0) + row.get("Pct_Aqua", 0)
        is_polder_zone = pct_polder > 50.0  # >50% rice/aqua = likely polder
        
        # Lateral flow for this segment
        Q_lat = base_flow_rate * area_km2
        
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
            "Is_Polder_Zone": is_polder_zone,
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


def add_point_sources(loads_df: pd.DataFrame, season: str = "dry") -> pd.DataFrame:
    """
    Add point source loads from major cities.
    
    Cities are treated as concentrated outfalls, not diffuse sources.
    Flow rate is independent of season (sewage is continuous).
    
    Args:
        loads_df: Existing loads DataFrame
        season: Season (affects background, not point sources)
    
    Returns:
        Updated DataFrame with point sources appended
    """
    point_records = []
    
    for city_name, city_data in POINT_SOURCES.items():
        branch = city_data["branch"]
        dist_km = city_data["distance_km"]
        pop = city_data["population"]
        
        segment_idx = int(dist_km * 1000 / DX_M)
        
        # Point source flow (sewage): ~200 L/person/day = 2.3e-6 m³/s/person
        Q_lat = pop * 200.0 / 86400.0 / 1000.0  # m³/s
        
        # Per-capita loads [g/s]
        NH4_load = pop * city_data["per_capita_NH4_g_day"] / 86400.0
        TOC_load = pop * city_data["per_capita_TOC_g_day"] / 86400.0
        NO3_load = NH4_load * 0.1  # Some nitrified
        PO4_load = NH4_load * 0.15  # N:P ~ 7:1
        DIC_load = TOC_load * 0.5   # Some respired
        
        point_records.append({
            "Branch": branch,
            "Segment_Index": segment_idx,
            "Distance_km": dist_km,
            "Area_km2": 0.0,  # Point source
            "Is_Polder_Zone": False,
            "Is_Point_Source": True,
            "Point_Source_Name": city_name,
            "Q_lat_m3_s": Q_lat,
            "NH4_load_g_s": NH4_load,
            "NO3_load_g_s": NO3_load,
            "PO4_load_g_s": PO4_load,
            "TOC_load_g_s": TOC_load,
            "DIC_load_g_s": DIC_load,
            "NH4_conc_mg_L": NH4_load / Q_lat if Q_lat > 0 else 0,
            "NO3_conc_mg_L": NO3_load / Q_lat if Q_lat > 0 else 0,
            "PO4_conc_mg_L": PO4_load / Q_lat if Q_lat > 0 else 0,
            "TOC_conc_mg_L": TOC_load / Q_lat if Q_lat > 0 else 0,
            "DIC_conc_mg_L": DIC_load / Q_lat if Q_lat > 0 else 0,
        })
    
    # Append point sources to main dataframe
    if point_records:
        point_df = pd.DataFrame(point_records)
        # Add Is_Point_Source column to original if missing
        if "Is_Point_Source" not in loads_df.columns:
            loads_df["Is_Point_Source"] = False
            loads_df["Point_Source_Name"] = ""
        return pd.concat([loads_df, point_df], ignore_index=True)
    
    return loads_df


def aggregate_to_grid(loads_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate 1-km segments to model grid cells (2 km spacing).
    
    The C-code needs loads at each grid cell, so we sum adjacent segments.
    Point sources are kept separate (not aggregated).
    """
    # Separate point sources
    is_point = loads_df.get("Is_Point_Source", pd.Series([False]*len(loads_df)))
    diffuse_df = loads_df[~is_point].copy()
    point_df = loads_df[is_point].copy()
    
    # Group diffuse loads by branch and segment index, sum the loads
    agg_df = diffuse_df.groupby(["Branch", "Segment_Index"]).agg({
        "Distance_km": "mean",
        "Area_km2": "sum",
        "Q_lat_m3_s": "sum",
        "NH4_load_g_s": "sum",
        "NO3_load_g_s": "sum",
        "PO4_load_g_s": "sum",
        "TOC_load_g_s": "sum",
        "DIC_load_g_s": "sum",
        "Is_Polder_Zone": "any",  # True if any segment is polder
    }).reset_index()
    
    # Recalculate concentrations for aggregated cells
    for sp in ["NH4", "NO3", "PO4", "TOC", "DIC"]:
        agg_df[f"{sp}_conc_mg_L"] = np.where(
            agg_df["Q_lat_m3_s"] > 1e-10,
            agg_df[f"{sp}_load_g_s"] / agg_df["Q_lat_m3_s"],
            0.0
        )
    
    # Add point sources back
    if len(point_df) > 0:
        # Ensure columns match
        for col in agg_df.columns:
            if col not in point_df.columns:
                point_df[col] = False if col == "Is_Polder_Zone" else 0.0
        agg_df = pd.concat([agg_df, point_df[agg_df.columns]], ignore_index=True)
    
    return agg_df


def generate_time_varying_loads(
    static_loads: pd.DataFrame,
    duration_days: float,
    output_path: Path
) -> None:
    """
    Generate time-varying lateral load files with Virtual Polder logic.
    
    This creates hourly load files that mimic gate operation:
    - Polder zones get flow=0 at high tide, flow=HIGH at low tide
    - Non-polder zones get constant flow
    - Point sources are always on (sewage is continuous)
    
    Args:
        static_loads: Base load DataFrame (from aggregate_to_grid)
        duration_days: Simulation duration [days]
        output_path: Directory to save time-varying files
    """
    print(f"\nGenerating time-varying loads with Virtual Polder logic...")
    print(f"  Duration: {duration_days} days")
    print(f"  Time step: {OUTPUT_TIME_STEP/3600:.1f} hours")
    
    # Generate tidal signal
    water_level = generate_tidal_signal(duration_days)
    gate_factor = calculate_polder_gate_factor(water_level)
    n_steps = len(water_level)
    
    # Time array in seconds
    time_s = np.arange(n_steps) * OUTPUT_TIME_STEP
    
    # Create output directory
    polder_dir = output_path / "polder_loads"
    polder_dir.mkdir(exist_ok=True)
    
    # Save time-varying factor file (for C-code to read)
    factor_df = pd.DataFrame({
        "Time_s": time_s,
        "Water_Level_m": water_level,
        "Gate_Factor": gate_factor,
    })
    factor_path = polder_dir / "gate_factor_timeseries.csv"
    factor_df.to_csv(factor_path, index=False)
    print(f"  Saved gate factor timeseries: {factor_path}")
    
    # Save polder zone identification (for C-code)
    polder_zones = static_loads[static_loads.get("Is_Polder_Zone", False) == True].copy()
    if len(polder_zones) > 0:
        polder_id_path = polder_dir / "polder_zones.csv"
        polder_zones[["Branch", "Segment_Index", "Distance_km"]].to_csv(
            polder_id_path, index=False
        )
        print(f"  Identified {len(polder_zones)} polder zone segments")
        print(f"  Saved polder zone list: {polder_id_path}")
    
    # Generate example time-varying load file for one branch
    # (Full implementation would do this for all branches)
    branches = static_loads["Branch"].unique()
    
    for branch in branches[:3]:  # Demo for first 3 branches
        branch_loads = static_loads[static_loads["Branch"] == branch].copy()
        
        # Create time-varying DataFrame
        records = []
        for t_idx, (t, gf) in enumerate(zip(time_s[:24*7], gate_factor[:24*7])):  # First week
            for _, row in branch_loads.iterrows():
                is_polder = row.get("Is_Polder_Zone", False)
                is_point = row.get("Is_Point_Source", False)
                
                if is_point:
                    # Point sources: constant
                    q_factor = 1.0
                elif is_polder:
                    # Polder zones: gate-controlled
                    q_factor = gf
                else:
                    # Non-polder diffuse: constant
                    q_factor = 1.0
                
                records.append({
                    "Time_s": t,
                    "Segment_Index": row["Segment_Index"],
                    "Q_lat_m3_s": row["Q_lat_m3_s"] * q_factor,
                    "NH4_load_g_s": row["NH4_load_g_s"] * q_factor,
                    "TOC_load_g_s": row["TOC_load_g_s"] * q_factor,
                })
        
        if records:
            tv_df = pd.DataFrame(records)
            tv_path = polder_dir / f"{branch}_timevarying_loads.csv"
            tv_df.to_csv(tv_path, index=False)
    
    print(f"  Generated time-varying load files for {min(3, len(branches))} branches")
    print(f"  Gate factor range: {gate_factor.min():.2f} - {gate_factor.max():.2f}")


def generate_full_year_loads(
    landuse_df: pd.DataFrame,
    output_path: Path,
    enable_polders: bool = True
) -> pd.DataFrame:
    """
    Generate a complete year of time-varying lateral loads with seasonal variation.
    
    This is the RECOMMENDED function for production runs, as it captures:
    - Full seasonal cycle (dry → wet → dry)
    - Monthly runoff variation based on MRC data
    - Virtual Polder gate operation
    - Spring-neap tidal modulation
    
    Output Files:
    - lateral_loads_annual.csv: Daily load factors for each month
    - lateral_loads_monthly/: Per-month static load files
    - polder_loads/: Hourly polder gate factors
    
    Args:
        landuse_df: Land use DataFrame
        output_path: Case directory path
        enable_polders: Enable Virtual Polder logic
    
    Returns:
        Summary DataFrame with annual totals
    """
    print("\n" + "=" * 70)
    print("GENERATING FULL YEAR LATERAL LOADS (365 DAYS)")
    print("=" * 70)
    
    # Create directories
    monthly_dir = output_path / "lateral_loads_monthly"
    monthly_dir.mkdir(exist_ok=True)
    
    annual_records = []
    monthly_summaries = []
    
    for month in range(1, 13):
        month_name = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][month - 1]
        
        # Determine season
        if month in DRY_MONTHS:
            season = "dry"
        elif month in [5, 11]:
            season = "transition"
        else:
            season = "wet"
        
        # Get monthly flow multiplier
        flow_mult = MONTHLY_FLOW_MULTIPLIERS[month]
        
        print(f"\n  Month {month:2d} ({month_name}): {season.upper():10s} | Flow mult: {flow_mult:.1f}x")
        
        # Calculate loads with seasonal flow rate
        base_flow_rate = get_seasonal_flow_rate("dry")  # Use dry as base
        monthly_flow_rate = base_flow_rate * flow_mult
        
        # Temporarily override the flow rate
        records = []
        for _, row in landuse_df.iterrows():
            branch = row["Branch"]
            dist_km = row["Distance_km"]
            area_km2 = row["Segment_Area_km2"]
            segment_idx = int(dist_km * 1000 / DX_M)
            
            pct_polder = row.get("Pct_Rice", 0) + row.get("Pct_Aqua", 0)
            is_polder_zone = pct_polder > 50.0
            
            Q_lat = monthly_flow_rate * area_km2
            
            # Calculate mass loads
            loads = {sp: 0.0 for sp in ["NH4", "NO3", "PO4", "TOC", "DIC"]}
            for lu_type in ["Urban", "Rice", "Aqua", "Mangrove", "Fruit"]:
                pct = row.get(f"Pct_{lu_type}", 0.0) / 100.0
                area_lu = area_km2 * pct
                for species, ef in EMISSION_FACTORS[lu_type].items():
                    mass_flux = ef * KG_HA_YR_TO_G_S_KM2 * area_lu
                    loads[species] += mass_flux
            
            # Scale loads by flow multiplier
            for sp in loads:
                loads[sp] *= flow_mult
            
            records.append({
                "Branch": branch,
                "Segment_Index": segment_idx,
                "Distance_km": dist_km,
                "Month": month,
                "Season": season,
                "Flow_Multiplier": flow_mult,
                "Q_lat_m3_s": Q_lat,
                "NH4_load_g_s": loads["NH4"],
                "NO3_load_g_s": loads["NO3"],
                "PO4_load_g_s": loads["PO4"],
                "TOC_load_g_s": loads["TOC"],
                "DIC_load_g_s": loads["DIC"],
                "Is_Polder_Zone": is_polder_zone,
            })
        
        month_df = pd.DataFrame(records)
        
        # Save monthly file
        month_path = monthly_dir / f"lateral_sources_month{month:02d}.csv"
        month_df.to_csv(month_path, index=False)
        
        # Compute monthly summary
        monthly_summaries.append({
            "Month": month,
            "Month_Name": month_name,
            "Season": season,
            "Flow_Multiplier": flow_mult,
            "Total_Q_m3_day": month_df["Q_lat_m3_s"].sum() * 86400,
            "Total_NH4_kg_day": month_df["NH4_load_g_s"].sum() * 86400 / 1000,
            "Total_NO3_kg_day": month_df["NO3_load_g_s"].sum() * 86400 / 1000,
            "Total_TOC_kg_day": month_df["TOC_load_g_s"].sum() * 86400 / 1000,
        })
        
        annual_records.extend(records)
    
    # Save annual summary
    summary_df = pd.DataFrame(monthly_summaries)
    summary_path = output_path / "lateral_loads_annual_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"\n  Saved annual summary: {summary_path}")
    
    # Generate daily interpolated file for C-code
    print("\n  Generating daily interpolated load factors...")
    daily_records = []
    for day in range(365):
        # Approximate month (0-indexed)
        month = (day // 30) + 1
        if month > 12:
            month = 12
        
        flow_mult = MONTHLY_FLOW_MULTIPLIERS[month]
        
        daily_records.append({
            "Day": day,
            "Month": month,
            "Flow_Multiplier": flow_mult,
        })
    
    daily_df = pd.DataFrame(daily_records)
    daily_path = output_path / "lateral_load_daily_factors.csv"
    daily_df.to_csv(daily_path, index=False)
    print(f"  Saved daily factors: {daily_path}")
    
    # Generate polder loads for full year if enabled
    if enable_polders:
        print("\n  Generating full-year polder gate factors...")
        generate_time_varying_loads(
            pd.DataFrame(annual_records[:len(landuse_df)]),  # Use first month as base
            FULL_YEAR_DAYS,
            output_path
        )
    
    # Print summary
    print("\n" + "=" * 70)
    print("ANNUAL LOAD SUMMARY")
    print("=" * 70)
    print(summary_df.to_string(index=False))
    
    # Annual totals
    print("\n" + "-" * 70)
    print("ANNUAL TOTALS:")
    print(f"  Total lateral flow:  {summary_df['Total_Q_m3_day'].mean() * 365:.0f} m³/year")
    print(f"  Total NH4 load:      {summary_df['Total_NH4_kg_day'].sum() * 30:.0f} kg N/year")
    print(f"  Total NO3 load:      {summary_df['Total_NO3_kg_day'].sum() * 30:.0f} kg N/year")
    print(f"  Total TOC load:      {summary_df['Total_TOC_kg_day'].sum() * 30:.0f} kg C/year")
    
    return summary_df


def main():
    """Main entry point with command-line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Generate lateral loads from land use map",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate full year data (RECOMMENDED for production)
  python generate_lateral_loads.py --annual
  
  # Dry season only
  python generate_lateral_loads.py --season dry
  
  # Wet season with 10x higher runoff
  python generate_lateral_loads.py --season wet
  
  # Enable Virtual Polder time-varying loads (30 days)
  python generate_lateral_loads.py --polders --duration 30
  
  # Custom case directory
  python generate_lateral_loads.py --case-dir INPUT/Cases/MyCase
        """
    )
    
    parser.add_argument(
        "--annual", "-a",
        action="store_true",
        help="Generate full year (365 days) of data covering all seasons (RECOMMENDED)"
    )
    parser.add_argument(
        "--season", "-s",
        choices=["dry", "wet", "transition"],
        default="dry",
        help="Season for runoff calculation (default: dry)"
    )
    parser.add_argument(
        "--polders", "-p",
        action="store_true",
        help="Enable Virtual Polder time-varying loads"
    )
    parser.add_argument(
        "--duration", "-d",
        type=float,
        default=365.0,
        help="Simulation duration in days for polder time series (default: 365)"
    )
    parser.add_argument(
        "--case-dir", "-c",
        type=str,
        default=None,
        help="Path to case directory (default: Mekong_Delta_Full)"
    )
    parser.add_argument(
        "--no-point-sources",
        action="store_true",
        help="Disable point source additions (cities)"
    )
    
    args = parser.parse_args()
    
    # Update case directory if specified
    global CASE_DIR
    if args.case_dir:
        CASE_DIR = Path(args.case_dir)
    
    print("=" * 70)
    print("GENERATING LATERAL LOADS FROM LAND USE MAP")
    print("=" * 70)
    
    # Load land use map
    print("\nLoading land use map...")
    landuse_df = load_landuse_map()
    print(f"  Loaded {len(landuse_df)} segments from {len(landuse_df['Branch'].unique())} branches")
    
    # ANNUAL MODE: Generate full year of data
    if args.annual:
        print("\n*** ANNUAL MODE: Generating full year (365 days) of data ***")
        summary = generate_full_year_loads(
            landuse_df, 
            CASE_DIR, 
            enable_polders=args.polders
        )
        
        # Also generate default lateral_sources.csv using average conditions
        print("\n  Generating default lateral_sources.csv (annual average)...")
        loads_df = calculate_loads(landuse_df, season="transition", enable_polders=args.polders)
        if not args.no_point_sources:
            loads_df = add_point_sources(loads_df, season="transition")
        grid_loads = aggregate_to_grid(loads_df)
        default_path = CASE_DIR / "lateral_sources.csv"
        grid_loads.to_csv(default_path, index=False)
        print(f"  Saved default: {default_path}")
        
        return grid_loads
    
    # SINGLE SEASON MODE
    print(f"  Season: {args.season.upper()}")
    print(f"  Flow rate: {get_seasonal_flow_rate(args.season):.4f} m³/s/km²")
    print(f"  Virtual Polders: {'ENABLED' if args.polders else 'disabled'}")
    
    # Calculate loads
    print("\nCalculating lateral loads...")
    loads_df = calculate_loads(landuse_df, season=args.season, enable_polders=args.polders)
    
    # Add point sources
    if not args.no_point_sources:
        print("\nAdding point sources (cities)...")
        loads_df = add_point_sources(loads_df, season=args.season)
        n_points = loads_df.get("Is_Point_Source", pd.Series([False]*len(loads_df))).sum()
        print(f"  Added {n_points} point sources")
    
    # Aggregate to grid
    print("\nAggregating to model grid (2 km cells)...")
    grid_loads = aggregate_to_grid(loads_df)
    
    # Save static output
    output_path = CASE_DIR / f"lateral_sources_{args.season}.csv"
    grid_loads.to_csv(output_path, index=False)
    print(f"\nSaved lateral sources to: {output_path}")
    
    # Also save as default name for backward compatibility
    default_path = CASE_DIR / "lateral_sources.csv"
    grid_loads.to_csv(default_path, index=False)
    
    # Generate time-varying loads if polders enabled
    if args.polders:
        generate_time_varying_loads(grid_loads, args.duration, CASE_DIR)
    
    # Print summary
    print("\n" + "=" * 70)
    print("LATERAL LOAD SUMMARY BY BRANCH")
    print("=" * 70)
    
    # Exclude point sources from diffuse summary
    is_point = grid_loads.get("Is_Point_Source", pd.Series([False]*len(grid_loads)))
    diffuse_loads = grid_loads[~is_point]
    
    summary = diffuse_loads.groupby("Branch").agg({
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
    
    print("\nDiffuse loads by branch (daily):")
    print(summary_daily[["Q_lat_m3_day", "NH4_kg_day", "NO3_kg_day", 
                         "PO4_kg_day", "TOC_kg_day", "DIC_kg_day"]].round(1).to_string())
    
    # Point source summary
    point_loads = grid_loads[is_point]
    if len(point_loads) > 0:
        print("\n" + "=" * 70)
        print("POINT SOURCE SUMMARY")
        print("=" * 70)
        for _, row in point_loads.iterrows():
            name = row.get("Point_Source_Name", "Unknown")
            print(f"  {name}:")
            print(f"    Flow: {row['Q_lat_m3_s']*86400:.0f} m³/day")
            print(f"    NH4:  {row['NH4_load_g_s']*86400/1000:.1f} kg/day")
            print(f"    TOC:  {row['TOC_load_g_s']*86400/1000:.1f} kg/day")
    
    # Identify polder zones
    polder_zones = grid_loads[grid_loads.get("Is_Polder_Zone", False) == True]
    if len(polder_zones) > 0:
        print("\n" + "=" * 70)
        print(f"POLDER ZONES ({len(polder_zones)} segments)")
        print("=" * 70)
        print("  These zones will have time-varying loads if --polders is enabled")
        for branch in polder_zones["Branch"].unique():
            n = len(polder_zones[polder_zones["Branch"] == branch])
            print(f"    {branch}: {n} segments")
    
    # Total delta loads
    print("\n" + "=" * 70)
    print(f"TOTAL MEKONG DELTA LATERAL LOADS ({args.season.upper()} SEASON)")
    print("=" * 70)
    
    totals = grid_loads[["Q_lat_m3_s", "NH4_load_g_s", "NO3_load_g_s", 
                         "PO4_load_g_s", "TOC_load_g_s", "DIC_load_g_s"]].sum()
    print(f"  Total lateral flow:  {totals['Q_lat_m3_s']*86400:.0f} m³/day")
    print(f"  Total NH4 load:      {totals['NH4_load_g_s']*86400/1000:.1f} kg N/day")
    print(f"  Total NO3 load:      {totals['NO3_load_g_s']*86400/1000:.1f} kg N/day")
    print(f"  Total PO4 load:      {totals['PO4_load_g_s']*86400/1000:.1f} kg P/day")
    print(f"  Total TOC load:      {totals['TOC_load_g_s']*86400/1000:.1f} kg C/day")
    print(f"  Total DIC load:      {totals['DIC_load_g_s']*86400/1000:.1f} kg C/day")
    
    return grid_loads


if __name__ == "__main__":
    main()
