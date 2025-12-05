#!/usr/bin/env python3
"""
Setup Mekong Delta - HIERARCHICAL Nguyen-Savenije Configuration
================================================================

This topology follows the approach of Nguyen et al. (2008) and Savenije (2005/2012):
- Merged distributaries for proper Savenije convergent geometry
- CORRECT hierarchical split: Co Chien branches first, then My Tho/Ham Luong
- 3 Tien outlets merged: Cua Tieu + Cua Dai → My Tho (North)
- 2 Co Chien merged: Cung Hau + Co Chien → Co Chien (South)  
- 2 Hau outlets merged: Dinh An + Tran De → Hau (single branch)

GEOGRAPHIC ORDER (North to South):
==================================
1. My Tho (formerly Cua Tieu/Cua Dai) - Northernmost Tien outlet
2. Ham Luong - Central Tien outlet
3. Co Chien (incl. Cung Hau) - Southernmost Tien outlet
4. Hau River (Dinh An + Tran De) - Single Hau outlet

CORRECT HIERARCHICAL TOPOLOGY:
==============================

    [Tan Chau]                      [Chau Doc]
       (1)                             (2)
        │                               │
    Tien_Main (50 km)              Hau_Main (45 km)
        │                               │
    VamNao_Tien (3)═══ Vam_Nao ═══VamNao_Hau (4)
        │               7 km            │
        │             PRISMATIC         │
    Tien_Lower (50 km)             Hau_River (90 km)
        │                               │
  Co_Chien_Split (5)              Hau_Mouth (9)
       / \\                             Sea
      /   \\
Co_Chien  Tien_Connector
  85km       15km
    │          │
CoChien_   MyTho_Split (10)
 Mouth        / \\
  (8)        /   \\
            /     \\
        My_Tho  Ham_Luong
         55km    60km
           │       │
      MyTho_   HamLuong_
       Mouth     Mouth
        (6)       (7)

SCIENTIFIC NOTES:
- Vam Nao: PRISMATIC (600m constant width) - pure hydraulic connector
- Each distributary: Convergent funnel shape (Savenije geometry)
- Node 5: Co_Chien_Split - Co Chien branches off first (geographically correct)
- Node 10: MyTho_Split - My Tho and Ham Luong separate further downstream

References:
- Nguyen, A.D., & Savenije, H.H.G. (2006). Salt intrusion in multi-channel estuaries
- Savenije, H.H.G. (2005, 2012). Salinity and Tides in Alluvial Estuaries

Author: Nguyen Truong An
Date: November 2024 (Updated: Hierarchical topology fix)
"""

import os
from pathlib import Path
from datetime import datetime

# ===========================================================================
# CONFIGURATION
# ===========================================================================

CASE_NAME = "Mekong_Delta_Full"
SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent
CASE_DIR = PROJECT_ROOT / "INPUT" / "Cases" / CASE_NAME
FORCING_DIR = CASE_DIR / "forcing_data"

# ===========================================================================
# HIERARCHICAL TOPOLOGY (10 Nodes, 9 Branches)
# ===========================================================================

# Node definitions:
# 1: Tan Chau (Tien inlet) - DISCHARGE BC
# 2: Chau Doc (Hau inlet) - DISCHARGE BC  
# 3: Vam Nao junction - Tien side
# 4: Vam Nao junction - Hau side (also Hau main stem continues)
# 5: Co_Chien_Split (first Tien distributary split - Co Chien branches off)
# 6: My Tho mouth (North) - LEVEL BC
# 7: Ham Luong mouth (Central) - LEVEL BC
# 8: Co Chien mouth (South) - LEVEL BC
# 9: Hau River mouth (Far South) - LEVEL BC
# 10: MyTho_Split (second Tien split - My Tho and Ham Luong separate)

NODES = {
    1: {"name": "Tan_Chau", "type": "DISCHARGE", "desc": "Tien inlet (60% Mekong)"},
    2: {"name": "Chau_Doc", "type": "DISCHARGE", "desc": "Hau inlet (40% Mekong)"},
    3: {"name": "VamNao_Tien", "type": "JUNCTION", "desc": "Vam Nao - Tien side"},
    4: {"name": "VamNao_Hau", "type": "JUNCTION", "desc": "Vam Nao - Hau side"},
    5: {"name": "Co_Chien_Split", "type": "JUNCTION", "desc": "First Tien split - Co Chien branches off"},
    6: {"name": "MyTho_Mouth", "type": "LEVEL", "desc": "My Tho (Cua Tieu/Dai) - Northernmost"},
    7: {"name": "HamLuong_Mouth", "type": "LEVEL", "desc": "Ham Luong - Central"},
    8: {"name": "CoChien_Mouth", "type": "LEVEL", "desc": "Co Chien (incl Cung Hau) - South"},
    9: {"name": "Hau_Mouth", "type": "LEVEL", "desc": "Hau River (Dinh An + Tran De) - Far South"},
    10: {"name": "MyTho_Split", "type": "JUNCTION", "desc": "Second Tien split - My Tho/Ham Luong separate"},
}

# Branch definitions following Nguyen-Savenije geometry
# Format: (ID, Name, NodeUp, NodeDown, Length_m, W_Up_m, W_Down_m, Depth_m, Chezy, Group, RS)
#
# KEY PHYSICS POINTS:
# 1. Vam_Nao: PRISMATIC (W_up = W_down) - No Savenije convergence
# 2. Distributaries: Strong convergence (W_down >> W_up) - Savenije funnel
# 3. Depths based on Nguyen et al. (2006) bathymetry
#
# HIERARCHICAL TOPOLOGY:
# - Co Chien branches off FIRST at Co_Chien_Split (Node 5)
# - My Tho and Ham Luong branch SECOND at MyTho_Split (Node 10)
# - Tien_Connector links the two split junctions (15 km)

BRANCHES = [
    # =========================================================================
    # MAIN STEMS (Group 1: River characteristics)
    # =========================================================================
    
    # 1. Tien Main: Tan Chau → Vam Nao
    #    50 km from border to Vam Nao junction
    (1, "Tien_Main", 1, 3, 50000, 1000, 1200, 15.0, 65, 1, 1.0),
    
    # 2. Hau Main: Chau Doc → Vam Nao Hau side
    #    45 km from border to Vam Nao junction
    (2, "Hau_Main", 2, 4, 45000, 900, 1100, 14.0, 65, 1, 1.0),
    
    # =========================================================================
    # VAM NAO CONNECTOR (Group 3: PRISMATIC - Critical physics!)
    # =========================================================================
    
    # 3. Vam Nao: Connects Tien (Node 3) to Hau (Node 4)
    #    ~7 km natural channel, very deep (18m), narrow (600m constant)
    #    PRISMATIC: W_up = W_down = 600m (No Savenije convergence!)
    (3, "Vam_Nao", 3, 4, 7000, 600, 600, 18.0, 75, 3, 1.0),
    
    # =========================================================================
    # TIEN LOWER: Vam Nao → Co_Chien_Split (Group 1)
    # =========================================================================
    
    # 4. Tien Lower: Continues from Vam Nao junction to first distributary split
    #    50 km reach before Co Chien branches off
    (4, "Tien_Lower", 3, 5, 50000, 1200, 1400, 14.0, 65, 1, 1.0),
    
    # =========================================================================
    # TIEN CONNECTOR: Co_Chien_Split → MyTho_Split (Group 1)
    # =========================================================================
    
    # 5. Tien_Connector: Short reach connecting the two split junctions
    #    15 km from Co_Chien_Split to MyTho_Split
    (5, "Tien_Connector", 5, 10, 15000, 1100, 1200, 13.0, 65, 1, 1.0),
    
    # =========================================================================
    # TIEN DISTRIBUTARIES (Group 2: Savenije estuarine funnel shapes)
    # Geographic order: North (My Tho) → Central (Ham Luong) → South (Co Chien)
    # =========================================================================
    
    # 6. Co_Chien: Co_Chien_Split (5) → Sea (8) - SOUTHERNMOST TIEN
    #    Merged: Cung Hau + Co Chien
    #    85 km, branches off FIRST, deepest Tien (13m)
    (6, "Co_Chien", 5, 8, 85000, 1100, 2800, 13.0, 60, 2, 1.2),
    
    # 7. My_Tho: MyTho_Split (10) → Sea (6) - NORTHERNMOST
    #    Merged: Cua Tieu + Cua Dai
    #    55 km (shorter now due to Tien_Connector), strong convergence
    (7, "My_Tho", 10, 6, 55000, 900, 2200, 11.0, 60, 2, 1.2),
    
    # 8. Ham_Luong: MyTho_Split (10) → Sea (7) - CENTRAL
    #    60 km, moderate convergence, deeper (12m)
    (8, "Ham_Luong", 10, 7, 60000, 1000, 2500, 12.0, 60, 2, 1.2),
    
    # =========================================================================
    # HAU DISTRIBUTARY (Group 2: Single merged branch)
    # =========================================================================
    
    # 9. Hau River: Vam Nao Hau (4) → Sea (9) - FAR SOUTH
    #    Merged: Dinh An + Tran De
    #    90 km long reach from Vam Nao to sea, deep (14m)
    (9, "Hau_River", 4, 9, 90000, 1200, 3000, 14.0, 58, 2, 1.3),
]

# Tidal phase progression (North to South, ~10° per ~25km coastline)
TIDE_FILES = {
    6: ("MyTho_Tide.csv", 0),       # My Tho: Reference phase
    7: ("HamLuong_Tide.csv", 8),    # Ham Luong: +8° (25 km south)
    8: ("CoChien_Tide.csv", 15),    # Co Chien: +15° (45 km south)
    9: ("Hau_Tide.csv", 30),        # Hau: +30° (80 km south)
}


# ===========================================================================
# FILE GENERATORS
# ===========================================================================

def generate_topology_csv() -> str:
    """Generate topology.csv with hierarchical Nguyen-Savenije structure."""
    lines = [
        "# ============================================================================",
        "# TOPOLOGY: Mekong Delta (Nguyen-Savenije Hierarchical)",
        "# ============================================================================",
        "# Following Nguyen et al. (2006) and Savenije (2005/2012) approach:",
        "# - Merged distributaries for proper Savenije convergent geometry",
        "# - My Tho = Cua Tieu + Cua Dai (North)",
        "# - Co Chien = Cung Hau + Co Chien (South Tien)",
        "# - Hau = Dinh An + Tran De (single Hau outlet)",
        "#",
        "# HIERARCHICAL TOPOLOGY:",
        "# - Co Chien branches off FIRST at Node 5 (Co_Chien_Split)",
        "# - Tien_Connector (15 km) links to Node 10 (MyTho_Split)",
        "# - My Tho and Ham Luong branch SECOND at Node 10",
        "#",
        "# CRITICAL PHYSICS:",
        "# - Vam Nao is PRISMATIC (W_up = W_down = 600m)",
        "# - This ensures hydraulic connector behavior, not estuarine",
        "#",
        "# Geographic order (N→S): My Tho → Ham Luong → Co Chien → Hau",
        "#",
        "# Format: ID, Name, NodeUp, NodeDown, Length, W_Up, W_Down, Depth, Chezy, Group, RS",
        "# ============================================================================",
        "",
    ]
    
    for b in BRANCHES:
        bid, name, node_up, node_down, length, w_up, w_down, depth, chezy, group, rs = b
        lines.append(f"{bid}, {name}, {node_up}, {node_down}, {length}, {w_up}, {w_down}, {depth}, {chezy}, {group}, {rs}")
    
    return "\n".join(lines)


def generate_boundary_map_csv() -> str:
    """Generate boundary_map.csv."""
    lines = [
        "# ============================================================================",
        "# BOUNDARY CONDITIONS: Mekong Delta (Simplified)",
        "# ============================================================================",
        "# 2 discharge inlets, 4 tidal outlets (N→S order)",
        "# Format: Node_ID, Type, FilePath [, SpeciesOption]",
        "# ============================================================================",
        "",
        "# UPSTREAM DISCHARGE (60/40 split)",
        "1, DISCHARGE, forcing_data/Tien_Inlet.csv",
        "1, SPECIES, species_river.csv, ALL",
        "",
        "2, DISCHARGE, forcing_data/Hau_Inlet.csv",
        "2, SPECIES, species_river.csv, ALL",
        "",
        "# DOWNSTREAM TIDAL (North to South)",
        "# My Tho - Northernmost",
        "6, LEVEL, forcing_data/MyTho_Tide.csv",
        "6, SPECIES, species_ocean.csv, ALL",
        "",
        "# Ham Luong - Central",
        "7, LEVEL, forcing_data/HamLuong_Tide.csv",
        "7, SPECIES, species_ocean.csv, ALL",
        "",
        "# Co Chien - South Tien",
        "8, LEVEL, forcing_data/CoChien_Tide.csv",
        "8, SPECIES, species_ocean.csv, ALL",
        "",
        "# Hau River - Far South (merged Dinh An + Tran De)",
        "9, LEVEL, forcing_data/Hau_Tide.csv",
        "9, SPECIES, species_ocean.csv, ALL",
    ]
    return "\n".join(lines)


def generate_case_config() -> str:
    """Generate case_config.txt with full year simulation."""
    return f"""# ============================================================================
# C-GEM Network: Mekong Delta (Nguyen-Savenije Hierarchical)
# ============================================================================
# Hierarchical topology following Nguyen et al. (2006):
# - 9 branches, 10 nodes
# - Merged distributaries for proper Savenije geometry
# - Prismatic Vam Nao connector
# - Co Chien splits first, then My Tho/Ham Luong
#
# FULL YEAR SIMULATION with seasonal variation:
# - Dry season (Dec-May): Low discharge, max salinity intrusion
# - Wet season (Jun-Nov): High discharge, minimal intrusion
#
# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
# ============================================================================

CaseName = {CASE_NAME}

# File paths
Topology = topology.csv
BoundaryMap = boundary_map.csv
BiogeoParams = biogeo_params.txt
OutputDir = OUTPUT/{CASE_NAME}

# Time parameters - FULL YEAR SIMULATION
StartDate = 2024-01-01
Duration = 365
Warmup = 30

# Numerical parameters
TimeStep = 60.0
DELXI = 2000.0

# Output
WriteNetCDF = 1
WriteCSV = 1
WriteReactionRates = 0

# Biogeochemistry mode
# simplified_mode = 1 enables the 80/20 parsimonious biogeochemistry
# suitable for data-sparse tropical systems
SimplifiedMode = 1
"""


def generate_biogeo_params() -> str:
    """Generate biogeo_params.txt with safe mode."""
    return """# ============================================================================
# C-GEM Biogeochemistry Parameters: Mekong Delta
# ============================================================================
# Safe mode configuration for stable simulation

# Environment
water_temp = 28.0
ws = 0.0005

# Light
I0 = 350.0
kd1 = 0.20
kd2_spm = 0.018
kd2_phy1 = 0.008
kd2_phy2 = 0.006

# Phytoplankton
alpha1 = 0.02
pbmax1 = 2.5
kexc1 = 0.10
kgrowth1 = 0.10
kmaint1 = 0.025
kmort1 = 0.08

alpha2 = 0.015
pbmax2 = 2.0
kexc2 = 0.10
kgrowth2 = 0.08
kmaint2 = 0.02
kmort2 = 0.06

# Nutrients
kdsi1 = 3.0
kn1 = 1.5
kpo41 = 0.3
kn2 = 2.0
kpo42 = 0.5

# Decomposition
kox = 0.15
kdenit = 0.04
knit = 0.12

# Half-saturation
ktox = 40.0
ko2 = 8.0
ko2_nit = 4.0
kno3 = 4.0
knh4 = 2.0
kino2 = 6.0

# Stoichiometry
redn = 0.151
redp = 0.0094
redsi = 0.12

# Gas exchange
pco2_atm = 420.0
wind_speed = 3.0
wind_coeff = 0.251
schmidt_exp = -0.5
current_k_factor = 0.30

# Benthic
benthic_resp_20C = 70.0
benthic_Q10 = 2.0
benthic_NH4_flux = 2.5
benthic_PO4_flux = 0.15

# GHG
N2O_yield_nit = 0.005
N2O_yield_denit = 0.015
benthic_CH4_flux = 150.0
benthic_N2O_flux = 8.0

# Safety flags
simplified_mode = 1
ghg_passive_mode = 1
skip_bacteria = 1
skip_multipool_oc = 1
skip_ghg_dynamics = 0
skip_p_adsorption = 1
"""


def generate_species_river() -> str:
    """Generate species_river.csv with time_s column covering full year."""
    # Generate a full year (365 days) of river species data
    # with seasonal variations in nutrients and phytoplankton
    
    lines = ["time_s,salinity,phy1,phy2,dsi,no3,nh4,po4,o2,toc,spm,dic,at,pco2,co2,ph,hs"]
    
    # Seasonal variations (based on MRC water quality data):
    # - Dry season (Jan-May): Lower discharge, higher concentrations
    # - Wet season (Jun-Nov): Higher discharge, diluted concentrations
    # - SPM peaks during flood season
    
    for day in range(0, 366, 1):  # Daily resolution for 365 days
        time_s = day * 86400
        month = (day // 30) + 1
        if month > 12:
            month = 12
        
        # Seasonal factors
        if month in [12, 1, 2, 3, 4, 5]:  # Dry season
            sal = 0.1
            phy1, phy2 = 20.0, 12.0  # Higher phyto in clear water
            dsi, no3, nh4, po4 = 100.0, 40.0, 15.0, 2.0  # Higher nutrients
            o2 = 240.0
            toc = 400.0
            spm = 50.0  # Low SPM
            dic, at = 1900.0, 2000.0
        elif month in [6, 7, 8, 9, 10]:  # Wet season
            sal = 0.1
            phy1, phy2 = 10.0, 5.0  # Lower phyto (turbidity limits light)
            dsi, no3, nh4, po4 = 70.0, 25.0, 8.0, 1.0  # Diluted nutrients
            o2 = 260.0
            toc = 250.0
            spm = 200.0 if month in [8, 9] else 120.0  # High SPM in flood
            dic, at = 1700.0, 1800.0
        else:  # Transition
            sal = 0.1
            phy1, phy2 = 15.0, 8.0
            dsi, no3, nh4, po4 = 85.0, 32.0, 11.0, 1.5
            o2 = 250.0
            toc = 320.0
            spm = 80.0
            dic, at = 1800.0, 1900.0
        
        lines.append(f"{time_s},{sal},{phy1},{phy2},{dsi},{no3},{nh4},{po4},"
                     f"{o2},{toc},{spm},{dic},{at},0.0,0.0,7.5,0.5")
    
    return "\n".join(lines)


def generate_species_ocean() -> str:
    """Generate species_ocean.csv with time_s column covering full year."""
    # Ocean conditions are more stable, but slight seasonal variations exist
    
    lines = ["time_s,salinity,phy1,phy2,dsi,no3,nh4,po4,o2,toc,spm,dic,at,pco2,co2,ph,hs"]
    
    for day in range(0, 366, 1):  # Daily resolution for 365 days
        time_s = day * 86400
        month = (day // 30) + 1
        if month > 12:
            month = 12
        
        # Slight seasonal variation in coastal waters
        if month in [12, 1, 2, 3, 4, 5]:  # Dry season - more marine influence
            sal = 33.0
            phy1, phy2 = 6.0, 4.0
            dsi, no3, nh4, po4 = 4.0, 1.5, 0.3, 0.2
            spm = 8.0
        else:  # Wet season - river plume influence
            sal = 28.0  # Lower salinity from river plume
            phy1, phy2 = 8.0, 5.0  # Higher coastal productivity
            dsi, no3, nh4, po4 = 8.0, 3.0, 0.8, 0.4
            spm = 20.0  # Higher from river input
        
        lines.append(f"{time_s},{sal},{phy1},{phy2},{dsi},{no3},{nh4},{po4},"
                     f"210.0,100.0,{spm},2100.0,2350.0,0.0,0.0,8.1,0.1")
    
    return "\n".join(lines)


def generate_hau_tide() -> str:
    """
    Generate Hau tide file (merged Dinh An + Tran De) for full year.
    
    PHASE 5 UPGRADE (December 2025 Audit):
    ======================================
    Added K1 (diurnal) and S2 (solar semidiurnal) constituents for journal-grade
    accuracy. The Mekong/South China Sea is a MIXED SEMIDIURNAL system where
    the diurnal inequality is significant.
    
    Tidal Constituents (from Nguyen et al. 2006, ADCP measurements):
    - M2: Principal lunar semidiurnal (12.42 h) - dominant
    - K1: Lunisolar diurnal (23.93 h) - creates high-high / low-low pattern
    - S2: Principal solar semidiurnal (12.00 h) - spring-neap modulation
    - O1: Lunar diurnal (25.82 h) - minor, combined with K1 here
    
    The form factor F = (K1 + O1) / (M2 + S2) ≈ 0.25-0.50 indicates mixed tide.
    
    Reference: Nguyen et al. (2006) Continental Shelf Research
    """
    import numpy as np
    
    # Full year at 30-min intervals
    dt = 1800  # seconds
    n_points = int(365 * 24 * 3600 / dt) + 1
    times = np.arange(n_points) * dt
    
    # ==========================================================================
    # TIDAL CONSTITUENT PARAMETERS (Mekong Delta - Hau River mouth)
    # ==========================================================================
    # Periods [seconds]
    M2_period = 12.42 * 3600    # 44712 s - Principal lunar semidiurnal
    K1_period = 23.93 * 3600    # 86148 s - Lunisolar diurnal
    S2_period = 12.00 * 3600    # 43200 s - Principal solar semidiurnal
    
    # Amplitudes [m] - from Nguyen et al. (2006) harmonic analysis at Dinh An
    M2_amp = 1.60       # Dominant semidiurnal (reduced from 2.1 to account for S2)
    K1_amp = 0.45       # Diurnal component (creates inequality)
    S2_amp = 0.50       # Solar semidiurnal (spring-neap with M2)
    
    # Phases [radians] - Greenwich phase converted to local
    M2_phase = 30 * np.pi / 180     # ~30° phase lag
    K1_phase = 120 * np.pi / 180    # Diurnal phase
    S2_phase = 50 * np.pi / 180     # Solar phase
    
    # Spring-neap modulation period (from M2-S2 beat frequency)
    # T_spring_neap = 1 / (1/M2_period - 1/S2_period) ≈ 14.77 days
    spring_neap_period = 14.77 * 86400
    
    # Seasonal modulation of mean sea level (monsoon setup)
    # Wet season: +0.1m setup (sea level rise from discharge), Dry season: -0.05m
    annual_period = 365.25 * 86400
    
    lines = ["time_s,water_level_m"]
    for i, t in enumerate(times):
        # =======================================================================
        # HARMONIC SUPERPOSITION: H(t) = Σ Aᵢ cos(ωᵢt + φᵢ)
        # =======================================================================
        
        # M2 component (principal lunar semidiurnal)
        H_M2 = M2_amp * np.cos(2 * np.pi * t / M2_period - M2_phase)
        
        # K1 component (diurnal - creates high-high / low-low inequality)
        H_K1 = K1_amp * np.cos(2 * np.pi * t / K1_period - K1_phase)
        
        # S2 component (solar semidiurnal - creates spring-neap with M2)
        H_S2 = S2_amp * np.cos(2 * np.pi * t / S2_period - S2_phase)
        
        # Seasonal mean sea level offset (peaks in September)
        seasonal_offset = 0.08 * np.sin(2 * np.pi * (t - 180*86400) / annual_period)
        
        # Total water level
        level = H_M2 + H_K1 + H_S2 + seasonal_offset
        
        lines.append(f"{int(t)},{level:.4f}")
    
    return "\n".join(lines)


def generate_full_year_discharge() -> tuple:
    """
    Generate full year discharge data for Tien and Hau rivers.
    
    Based on MRC data, the Mekong has strong seasonal variation:
    - Dry season (Jan-May): ~2,000 m³/s at Tan Chau
    - Wet season (Jul-Oct): ~20,000 m³/s peak
    - Annual average: ~13,000 m³/s
    
    Returns:
        Tuple of (tien_discharge_csv, hau_discharge_csv)
    """
    import numpy as np
    
    # Full year at hourly intervals
    dt = 3600  # seconds
    n_points = int(365 * 24 * 3600 / dt) + 1
    times = np.arange(n_points) * dt
    
    # Monthly discharge multipliers (relative to mean)
    # Based on MRC data at Tan Chau
    monthly_Q_factors = {
        1: 0.25,   # January - Dry
        2: 0.20,   # February - Driest
        3: 0.18,   # March - Driest
        4: 0.22,   # April - Late dry
        5: 0.40,   # May - Early transition
        6: 0.80,   # June - Early wet
        7: 1.20,   # July - Wet
        8: 1.60,   # August - Peak flood
        9: 1.80,   # September - Peak flood
        10: 1.50,  # October - Late wet
        11: 0.80,  # November - Transition
        12: 0.40,  # December - Early dry
    }
    
    # Annual mean discharge (total Mekong at Cambodian border)
    Q_mean_total = 13000.0  # m³/s
    
    # Split: Tien 60%, Hau 40%
    Q_mean_tien = Q_mean_total * 0.60
    Q_mean_hau = Q_mean_total * 0.40
    
    tien_lines = ["time_s,discharge_m3_s"]
    hau_lines = ["time_s,discharge_m3_s"]
    
    for i, t in enumerate(times):
        # Determine month (approximate)
        day = t / 86400
        month = int((day % 365) / 30.44) + 1
        if month > 12:
            month = 12
        
        factor = monthly_Q_factors[month]
        
        # Add small daily variation (tidal influence on discharge)
        tidal_var = 0.05 * np.sin(2 * np.pi * t / (12.42 * 3600))
        
        Q_tien = Q_mean_tien * factor * (1.0 + tidal_var)
        Q_hau = Q_mean_hau * factor * (1.0 + tidal_var)
        
        tien_lines.append(f"{int(t)},{Q_tien:.1f}")
        hau_lines.append(f"{int(t)},{Q_hau:.1f}")
    
    return "\n".join(tien_lines), "\n".join(hau_lines)


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    """Generate all configuration files for full-year simulation."""
    import numpy as np
    
    print("=" * 70)
    print("Setting up Mekong Delta (Nguyen-Savenije Simplified)")
    print("FULL YEAR (365 days) with Seasonal Variation")
    print("=" * 70)
    print(f"Case directory: {CASE_DIR}")
    
    # Create directories
    print("\n1. Creating directories...")
    CASE_DIR.mkdir(parents=True, exist_ok=True)
    FORCING_DIR.mkdir(parents=True, exist_ok=True)
    
    # Generate files
    print("\n2. Generating topology.csv...")
    with open(CASE_DIR / "topology.csv", "w") as f:
        f.write(generate_topology_csv())
    print(f"   ✓ {len(BRANCHES)} branches, 10 nodes")
    print("   ✓ Vam Nao: PRISMATIC (600m x 600m)")
    print("   ✓ Hierarchical: Co_Chien splits first, then My_Tho/Ham_Luong")
    
    print("\n3. Generating boundary_map.csv...")
    with open(CASE_DIR / "boundary_map.csv", "w") as f:
        f.write(generate_boundary_map_csv())
    
    print("\n4. Generating case_config.txt (365 days)...")
    with open(CASE_DIR / "case_config.txt", "w") as f:
        f.write(generate_case_config())
    
    print("\n5. Generating biogeo_params.txt...")
    with open(CASE_DIR / "biogeo_params.txt", "w") as f:
        f.write(generate_biogeo_params())
    
    print("\n6. Generating species files (365 days with seasonal variation)...")
    with open(CASE_DIR / "species_river.csv", "w") as f:
        f.write(generate_species_river())
    with open(CASE_DIR / "species_ocean.csv", "w") as f:
        f.write(generate_species_ocean())
    
    print("\n7. Generating tidal files (365 days)...")
    with open(FORCING_DIR / "Hau_Tide.csv", "w") as f:
        f.write(generate_hau_tide())
    
    # Generate additional tidal files with phase offsets
    M2_period = 12.42 * 3600
    dt = 1800
    n_points = int(365 * 24 * 3600 / dt) + 1
    times = np.arange(n_points) * dt
    spring_neap_period = 14.77 * 86400
    annual_period = 365.25 * 86400
    
    for node_id, (filename, phase_deg) in TIDE_FILES.items():
        if filename == "Hau_Tide.csv":
            continue  # Already generated
        
        phase = phase_deg * np.pi / 180
        amplitude = 2.0  # Slightly lower for other outlets
        
        lines = ["time_s,water_level_m"]
        for t in times:
            spring_neap = 1.0 + 0.2 * np.cos(2 * np.pi * t / spring_neap_period)
            seasonal_offset = 0.05 * np.sin(2 * np.pi * (t - 180*86400) / annual_period)
            level = amplitude * spring_neap * np.cos(2 * np.pi * t / M2_period - phase) + seasonal_offset
            lines.append(f"{int(t)},{level:.4f}")
        
        with open(FORCING_DIR / filename, "w") as f:
            f.write("\n".join(lines))
        print(f"   ✓ {filename} (phase: {phase_deg}°)")
    
    print("\n8. Generating discharge files (365 days with seasonal variation)...")
    tien_csv, hau_csv = generate_full_year_discharge()
    with open(FORCING_DIR / "Tien_Inlet.csv", "w") as f:
        f.write(tien_csv)
    with open(FORCING_DIR / "Hau_Inlet.csv", "w") as f:
        f.write(hau_csv)
    print("   ✓ Tien_Inlet.csv (60% of Mekong flow)")
    print("   ✓ Hau_Inlet.csv (40% of Mekong flow)")
    
    # Print topology diagram
    print("\n" + "=" * 70)
    print("HIERARCHICAL NETWORK TOPOLOGY (Nguyen-Savenije)")
    print("=" * 70)
    print("""
                    [Tan Chau]              [Chau Doc]
                       (1)                     (2)
                        │                       │
                    Tien_Main               Hau_Main
                      50 km                  45 km
                        │                       │
    [Vam Nao Jct] ════(3)════ Vam_Nao ════(4)════
       Tien side       │      7 km PRISMATIC    │
                       │                        │
                   Tien_Lower               Hau_River
                     50 km                    90 km
                       │                        │
              [Co_Chien_Split]                  │
                    (5)                         │
                   /   \\                        │
                  /     \\                       │
            Co_Chien  Tien_Connector            │
              85km       15km                   │
                │          │                    │
              (8)    [MyTho_Split]              │
             [Sea]        (10)                  │
                         /   \\                  │
                        /     \\                 │
                    My_Tho  Ham_Luong           │
                     55km    60km               │
                       │       │                │
                     (6)     (7)              (9)
                    [Sea]   [Sea]            [Sea]
                      ↑       ↑                ↑
                   NORTH   CENTER          FAR SOUTH

    Hierarchical split order:
    1. Co_Chien branches FIRST at Node 5 (Co_Chien_Split)
    2. Tien_Connector (15km) continues to Node 10 (MyTho_Split)  
    3. My_Tho and Ham_Luong branch SECOND at Node 10
    
    Geographic order: My Tho (N) → Ham Luong → Co Chien → Hau (S)
    """)
    
    print("=" * 70)
    print("SEASONAL VARIATION SUMMARY")
    print("=" * 70)
    print("""
    Month   Season       Q_total   Salinity Intrusion
    -----   ------       -------   ------------------
    Jan     Dry          2,600     Maximum (60-70 km)
    Feb     Dry (min)    2,340     Maximum
    Mar     Dry (min)    2,080     Maximum
    Apr     Dry          2,860     High
    May     Transition   5,200     Decreasing
    Jun     Early Wet    10,400    Moderate
    Jul     Wet          15,600    Low
    Aug     Peak Flood   20,800    Minimal (<20 km)
    Sep     Peak Flood   23,400    Minimal
    Oct     Late Wet     19,500    Low
    Nov     Transition   10,400    Increasing
    Dec     Early Dry    5,200     High
    """)
    
    print("=" * 70)
    print("Setup Complete!")
    print("=" * 70)
    print("\nNext steps:")
    print("  1. Generate land use: python scripts/generate_synthetic_landuse.py")
    print("  2. Generate lateral loads: python scripts/generate_lateral_loads.py --annual")
    print("  3. Build: scripts/build.bat")
    print("  4. Run: bin/Debug/CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt")
    print("\n  Expected runtime for 365-day simulation: ~10-15 minutes")
    
    return 0


if __name__ == "__main__":
    exit(main())
