#!/usr/bin/env python3
"""
Setup Mekong Delta HIERARCHICAL Case (Correct Physics)
======================================================

Generates input files for the ACADEMIC Mekong Delta case.
This script implements the CORRECT topology based on actual Mekong geography.

HIERARCHICAL TOPOLOGY (North to South):
========================================

TIEN RIVER (Northern Branch - 60% of Mekong flow):
  1. Tan Chau (Node 1) -> Vam Nao Junction (Node 3)
  2. Node 3 -> My Thuan (Node 5)
     -> SPLIT: Co Chien leaves here (southernmost Tien distributary)
  3. Node 5 -> Cho Lach (Node 6)
     -> SPLIT: Ham Luong leaves here (central Tien distributary)
  4. Node 6 -> My Tho mouth (Node 8) (northernmost Tien distributary)

HAU RIVER (Southern Branch - 40% of Mekong flow):
  1. Chau Doc (Node 2) -> Vam Nao Junction (Node 4)
  2. Node 4 -> Dai Ngai/Can Tho (Node 7)
     -> SPLIT: Dinh An and Tran De

VAM NAO CONNECTION (Critical for flow redistribution):
  Node 3 (Tien) <-> Node 4 (Hau)
  
  **PHYSICS KEY**: Defined as PRISMATIC (Width_Up == Width_Down)
  - This forces the Saint-Venant solver to treat it as a hydraulic shunt
  - NOT an estuary with Savenije convergence physics
  - Flow driven purely by head difference (dH) between Tien and Hau

OUTLET ORDER (North to South, matching tidal phase progression):
  Node 8:  My Tho mouth (Cua Tieu/Dai) - Phase 0°
  Node 9:  Ham Luong mouth            - Phase 10°
  Node 10: Co Chien mouth             - Phase 20°
  Node 11: Dinh An mouth              - Phase 30°
  Node 12: Tran De mouth              - Phase 40°

References:
- Nguyen et al. (2006) - Mekong Delta geometry
- Vo Quoc Thanh (2021) - Updated bathymetry
- Fujii et al. (2003) - Vam Nao hydraulics
- MRC Hydrological Database - Flow splits

Author: CGEM Development Team
Date: November 2024
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
# NODE DEFINITIONS (12 Nodes Total)
# ===========================================================================
# Nodes 1-2:   Upstream discharge boundaries
# Nodes 3-7:   Internal junctions
# Nodes 8-12:  Ocean tidal boundaries (North to South)

NODES = {
    # Upstream Inlets
    1: {"name": "Tan_Chau", "type": "DISCHARGE", "desc": "Tien River inlet (60% Mekong)"},
    2: {"name": "Chau_Doc", "type": "DISCHARGE", "desc": "Hau River inlet (40% Mekong)"},
    
    # Internal Junctions
    3: {"name": "VamNao_Tien", "type": "JUNCTION", "desc": "Vam Nao junction - Tien side"},
    4: {"name": "VamNao_Hau", "type": "JUNCTION", "desc": "Vam Nao junction - Hau side"},
    5: {"name": "My_Thuan", "type": "JUNCTION", "desc": "My Thuan - Co Chien splits off"},
    6: {"name": "Cho_Lach", "type": "JUNCTION", "desc": "Cho Lach - Ham Luong splits off"},
    7: {"name": "Dai_Ngai", "type": "JUNCTION", "desc": "Dai Ngai/Can Tho - Hau splits"},
    
    # Ocean Outlets (Ordered North to South)
    8:  {"name": "MyTho_Mouth", "type": "LEVEL", "desc": "My Tho (Cua Tieu) - Northernmost"},
    9:  {"name": "HamLuong_Mouth", "type": "LEVEL", "desc": "Ham Luong mouth"},
    10: {"name": "CoChien_Mouth", "type": "LEVEL", "desc": "Co Chien mouth"},
    11: {"name": "DinhAn_Mouth", "type": "LEVEL", "desc": "Dinh An - Main Hau mouth"},
    12: {"name": "TranDe_Mouth", "type": "LEVEL", "desc": "Tran De - Southernmost"},
}

# ===========================================================================
# BRANCH DEFINITIONS (11 Branches)
# ===========================================================================
# Format: (ID, Name, NodeUp, NodeDown, Length_m, W_Up_m, W_Down_m, Depth_m, Chezy, Group, RS)
#
# Group Classification:
#   1 = Main river channels (slight convergence)
#   2 = Estuarine distributaries (Savenije funnel shape)
#   3 = Vam Nao connector (PRISMATIC - critical for correct physics)
#
# Width Convention:
#   - Converging (W_Up < W_Down): Estuarine Savenije physics
#   - Prismatic (W_Up == W_Down): Standard hydraulic channel physics

BRANCHES = [
    # =========================================================================
    # UPSTREAM SEGMENTS (Main River Characteristics)
    # =========================================================================
    
    # 1. Tien Upper: Tan Chau (1) -> Vam Nao Junction (3)
    #    Main stem before diversion, relatively narrow
    (1, "Tien_Up", 1, 3, 50000, 1000, 1200, 15.0, 65, 1, 1.0),
    
    # 2. Hau Upper: Chau Doc (2) -> Vam Nao Junction (4)
    #    Parallel main stem, slightly narrower than Tien
    (2, "Hau_Up", 2, 4, 45000, 900, 1100, 14.0, 65, 1, 1.0),
    
    # =========================================================================
    # VAM NAO CONNECTION (THE KEY PHYSICS FIX!)
    # =========================================================================
    
    # 3. Vam Nao: Tien side (3) -> Hau side (4)
    #    **CRITICAL**: PRISMATIC geometry (600m -> 600m)
    #    - Deep channel (18m) allows gravitational exchange
    #    - High Chezy (75) = low friction for efficient water transfer
    #    - Prismatic shape means NO Savenije convergence terms
    #    - Flow driven purely by head difference dH(Tien) - dH(Hau)
    (3, "Vam_Nao", 3, 4, 7000, 600, 600, 18.0, 75, 3, 1.0),
    
    # =========================================================================
    # TIEN MIDDLE SEGMENTS (Hierarchical Splits)
    # =========================================================================
    
    # 4. Tien Mid A: Vam Nao Jct (3) -> My Thuan (5)
    #    After Vam Nao diversion, before first split
    (4, "Tien_Mid_A", 3, 5, 50000, 1200, 1400, 14.0, 65, 1, 1.0),
    
    # 5. Co Chien: My Thuan (5) -> Ocean (10)
    #    FIRST SPLIT from Tien - southernmost Tien distributary
    #    Strong funnel shape for Savenije estuarine physics
    (5, "Co_Chien", 5, 10, 85000, 1100, 2800, 13.0, 60, 2, 1.2),
    
    # 6. Tien Mid B: My Thuan (5) -> Cho Lach (6)
    #    Connecting segment between My Thuan and Cho Lach splits
    (6, "Tien_Mid_B", 5, 6, 25000, 1000, 1200, 13.0, 65, 1, 1.0),
    
    # 7. Ham Luong: Cho Lach (6) -> Ocean (9)
    #    SECOND SPLIT from Tien - central Tien distributary
    (7, "Ham_Luong", 6, 9, 75000, 1000, 2500, 12.0, 60, 2, 1.2),
    
    # 8. My Tho: Cho Lach (6) -> Ocean (8)
    #    FINAL Tien outlet - northernmost distributary (Cua Tieu/Dai)
    (8, "My_Tho", 6, 8, 70000, 900, 2200, 11.0, 60, 2, 1.2),
    
    # =========================================================================
    # HAU SEGMENTS (Simpler Split Structure)
    # =========================================================================
    
    # 9. Hau Middle: Vam Nao Jct (4) -> Dai Ngai (7)
    #    Main Hau stem after receiving Vam Nao inflow
    #    Longer reach before splitting
    (9, "Hau_Mid", 4, 7, 90000, 1100, 1500, 15.0, 65, 1, 1.0),
    
    # 10. Dinh An: Dai Ngai (7) -> Ocean (11)
    #     Main Hau distributary (largest mouth)
    (10, "Dinh_An", 7, 11, 45000, 1400, 3200, 12.0, 58, 2, 1.3),
    
    # 11. Tran De: Dai Ngai (7) -> Ocean (12)
    #     Secondary Hau distributary (southernmost outlet)
    (11, "Tran_De", 7, 12, 40000, 1100, 2500, 11.0, 58, 2, 1.3),
]

# ===========================================================================
# BOUNDARY CONFIGURATION
# ===========================================================================

DISCHARGE_NODES = [1, 2]
LEVEL_NODES = [8, 9, 10, 11, 12]  # North to South

# Forcing files (names match generate_synthetic_mekong.py output)
DISCHARGE_FILES = {
    1: "Tien_Inlet.csv",
    2: "Hau_Inlet.csv",
}

# Tidal files with phase progression (North to South)
TIDE_FILES = {
    8:  ("MyTho_Tide.csv", 0),      # Phase 0° - reference
    9:  ("HamLuong_Tide.csv", 10),  # Phase 10°
    10: ("CoChien_Tide.csv", 20),   # Phase 20°
    11: ("DinhAn_Tide.csv", 30),    # Phase 30°
    12: ("TranDe_Tide.csv", 40),    # Phase 40°
}


# ===========================================================================
# FILE GENERATORS
# ===========================================================================

def generate_topology_csv() -> str:
    """Generate topology.csv with correct hierarchical structure."""
    lines = [
        "# ============================================================================",
        "# TOPOLOGY: Mekong Delta Full Network (Hierarchical)",
        "# ============================================================================",
        "# Academically correct representation with:",
        "#   - Hierarchical splits (My Thuan -> Cho Lach)",
        "#   - PRISMATIC Vam Nao (correct hydraulic physics)",
        "#   - North-to-South outlet ordering",
        "#",
        "# Node Convention (1-12, contiguous):",
        "#   1-2:   Upstream discharge BCs (Tien, Hau)",
        "#   3-7:   Internal junctions",
        "#   8-12:  Ocean outlets (North to South)",
        "#",
        "# CRITICAL PHYSICS NOTE:",
        "#   Vam Nao is PRISMATIC (W_up = W_down = 600m)",
        "#   This ensures Saint-Venant solver treats it as hydraulic shunt,",
        "#   NOT an estuary with Savenije convergence.",
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
    """Generate boundary_map.csv with proper node mappings."""
    lines = [
        "# ============================================================================",
        "# BOUNDARY CONDITIONS: Mekong Delta Full (Hierarchical)",
        "# ============================================================================",
        "# Node mapping:",
        "#   1: Tan Chau (Tien inlet) - 60% Mekong discharge",
        "#   2: Chau Doc (Hau inlet)  - 40% Mekong discharge",
        "#   8-12: Ocean outlets (North to South)",
        "#",
        "# Format: Node_ID, Type, FilePath [, SpeciesOption]",
        "# ============================================================================",
        "",
        "# ===========================================================================",
        "# UPSTREAM DISCHARGE BOUNDARIES",
        "# ===========================================================================",
        "",
    ]
    
    # Discharge BCs
    for node_id in DISCHARGE_NODES:
        filename = DISCHARGE_FILES[node_id]
        node_info = NODES[node_id]
        lines.append(f"# Node {node_id}: {node_info['name']} - {node_info['desc']}")
        lines.append(f"{node_id}, DISCHARGE, forcing_data/{filename}")
        lines.append(f"{node_id}, SPECIES, species_river.csv, ALL")
        lines.append("")
    
    lines.append("# ===========================================================================")
    lines.append("# DOWNSTREAM TIDAL BOUNDARIES (North to South)")
    lines.append("# ===========================================================================")
    lines.append("")
    
    # Tidal BCs (ordered North to South)
    for node_id in LEVEL_NODES:
        filename, phase = TIDE_FILES[node_id]
        node_info = NODES[node_id]
        lines.append(f"# Node {node_id}: {node_info['name']} - Phase {phase}°")
        lines.append(f"{node_id}, LEVEL, forcing_data/{filename}")
        lines.append(f"{node_id}, SPECIES, species_ocean.csv, ALL")
        lines.append("")
    
    return "\n".join(lines)


def generate_case_config() -> str:
    """Generate case_config.txt with optimized numerical parameters."""
    return f"""# ============================================================================
# C-GEM Network Case Configuration: {CASE_NAME} (Hierarchical)
# ============================================================================
# Complete Mekong Delta with academically correct topology:
# - Hierarchical Tien splits (My Thuan -> Cho Lach)
# - PRISMATIC Vam Nao for correct hydraulic physics
# - 5 ocean outlets ordered North to South
#
# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
# ============================================================================

CaseName = {CASE_NAME}

# ---------------------------------------------------------------------------
# FILE PATHS
# ---------------------------------------------------------------------------
Topology = topology.csv
BoundaryMap = boundary_map.csv
BiogeoParams = biogeo_params.txt
OutputDir = OUTPUT/{CASE_NAME}

# ---------------------------------------------------------------------------
# TIME PARAMETERS
# ---------------------------------------------------------------------------
StartDate = 2024-03-01
Duration = 10
Warmup = 3

# ---------------------------------------------------------------------------
# NUMERICAL PARAMETERS
# ---------------------------------------------------------------------------
# Time step [seconds]
# Reduced to 60s for stability with the looped network topology
# CFL: dt < dx / (sqrt(g*h) + |u|) ~ 2000/(12+2) ~ 140s
TimeStep = 60.0

# Grid spacing [meters]
DELXI = 2000.0

# ---------------------------------------------------------------------------
# OUTPUT OPTIONS
# ---------------------------------------------------------------------------
WriteNetCDF = 1
WriteCSV = 1
WriteReactionRates = 0
"""


def generate_biogeo_params() -> str:
    """Generate biogeo_params.txt with safe mode settings."""
    return """# ============================================================================
# C-GEM Biogeochemistry Parameters: Mekong Delta (Safe Mode)
# ============================================================================
# Simplified configuration for stable simulation
# Uses 80/20 approach - bulk parameters for data-sparse regions
#
# SAFETY FLAGS:
#   simplified_mode = 1  -> Use bulk decay rates (not multi-pool)
#   ghg_passive_mode = 1 -> GHG as diagnostic only (stable)
# ============================================================================

# ---------------------------------------------------------------------------
# ENVIRONMENTAL CONDITIONS
# ---------------------------------------------------------------------------
water_temp = 28.0         # Mean temperature [°C] (tropical, year-round)
ws = 0.0005               # Settling velocity [m/s]

# ---------------------------------------------------------------------------
# LIGHT PARAMETERS
# ---------------------------------------------------------------------------
I0 = 350.0                # Surface irradiance [W/m²]
kd1 = 0.20                # Base light attenuation [1/m]
kd2_spm = 0.018           # SPM contribution to attenuation [m²/mg]
kd2_phy1 = 0.008          # PHY1 self-shading
kd2_phy2 = 0.006          # PHY2 self-shading

# ---------------------------------------------------------------------------
# PHYTOPLANKTON PARAMETERS
# ---------------------------------------------------------------------------
# PHY1 (Diatoms) - silica-limited in estuaries
alpha1 = 0.02             # P-I curve slope [µg C/µg Chl/W/m²]
pbmax1 = 2.5              # Max photosynthesis [1/day]
kexc1 = 0.10              # Excretion fraction
kgrowth1 = 0.10           # Growth respiration
kmaint1 = 0.025           # Maintenance [1/day]
kmort1 = 0.08             # Mortality [1/day]

# PHY2 (Green algae) - nitrogen-limited
alpha2 = 0.015
pbmax2 = 2.0
kexc2 = 0.10
kgrowth2 = 0.08
kmaint2 = 0.02
kmort2 = 0.06

# ---------------------------------------------------------------------------
# NUTRIENT HALF-SATURATION CONSTANTS [µmol/L]
# ---------------------------------------------------------------------------
kdsi1 = 3.0               # DSi for PHY1
kn1 = 1.5                 # N for PHY1
kpo41 = 0.3               # PO4 for PHY1
kn2 = 2.0                 # N for PHY2
kpo42 = 0.5               # PO4 for PHY2

# ---------------------------------------------------------------------------
# DECOMPOSITION RATES [1/day]
# ---------------------------------------------------------------------------
kox = 0.15                # Aerobic TOC mineralization
kdenit = 0.04             # Denitrification
knit = 0.12               # Nitrification (NH4 -> NO3)

# ---------------------------------------------------------------------------
# HALF-SATURATION CONSTANTS [µmol/L]
# ---------------------------------------------------------------------------
ktox = 40.0               # TOC for decomposition
ko2 = 8.0                 # O2 for aerobic respiration
ko2_nit = 4.0             # O2 for nitrification
kno3 = 4.0                # NO3 for denitrification
knh4 = 2.0                # NH4 for nitrification
kino2 = 6.0               # O2 inhibition threshold for denitrification

# ---------------------------------------------------------------------------
# STOICHIOMETRY (Redfield ratios)
# ---------------------------------------------------------------------------
redn = 0.151              # N:C molar ratio (16/106)
redp = 0.0094             # P:C molar ratio (1/106)
redsi = 0.12              # Si:C for diatoms

# ---------------------------------------------------------------------------
# GAS EXCHANGE
# ---------------------------------------------------------------------------
pco2_atm = 420.0          # Atmospheric pCO2 [µatm]
wind_speed = 3.0          # Mean wind [m/s]
wind_coeff = 0.251        # Wanninkhof coefficient
schmidt_exp = -0.5
current_k_factor = 0.30

# ---------------------------------------------------------------------------
# BENTHIC FLUXES
# ---------------------------------------------------------------------------
benthic_resp_20C = 70.0   # Benthic CO2 flux [mmol C/m²/day]
benthic_Q10 = 2.0
benthic_NH4_flux = 2.5    # Benthic NH4 flux [mmol N/m²/day]
benthic_PO4_flux = 0.15   # Benthic PO4 flux [mmol P/m²/day]

# ---------------------------------------------------------------------------
# GHG PARAMETERS (Passive diagnostic mode)
# ---------------------------------------------------------------------------
N2O_yield_nit = 0.005     # N2O yield from nitrification
N2O_yield_denit = 0.015   # N2O yield from denitrification
benthic_CH4_flux = 150.0  # Benthic CH4 flux [µmol/m²/day]
benthic_N2O_flux = 8.0    # Benthic N2O flux [nmol/m²/day]

# ============================================================================
# SAFETY MODE FLAGS (CRITICAL FOR STABILITY!)
# ============================================================================
simplified_mode = 1       # 1 = Use simplified bulk kinetics
ghg_passive_mode = 1      # 1 = GHG diagnostic only (no feedback)
skip_bacteria = 1         # 1 = Skip explicit bacteria dynamics
skip_multipool_oc = 1     # 1 = Use single TOC pool
skip_ghg_dynamics = 0     # 0 = Still calculate GHG (but passive)
skip_p_adsorption = 1     # 1 = Skip PO4-PIP equilibrium
"""


def generate_species_river() -> str:
    """Generate species_river.csv - upstream boundary conditions with time_s."""
    # Create time-series format that the loader expects
    return """time_s,salinity,phy1,phy2,dsi,no3,nh4,po4,o2,toc,spm,dic,at,pco2,co2,ph,hs
0,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
86400,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
172800,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
259200,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
345600,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
432000,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
518400,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
604800,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
691200,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
777600,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
864000,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
950400,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
1036800,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
1123200,0.1,15.0,8.0,80.0,35.0,12.0,1.5,250.0,300.0,80.0,1800.0,1900.0,0.0,0.0,7.5,0.5
"""


def generate_species_ocean() -> str:
    """Generate species_ocean.csv - downstream boundary conditions with time_s."""
    # Create time-series format that the loader expects
    return """time_s,salinity,phy1,phy2,dsi,no3,nh4,po4,o2,toc,spm,dic,at,pco2,co2,ph,hs
0,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
86400,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
172800,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
259200,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
345600,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
432000,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
518400,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
604800,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
691200,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
777600,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
864000,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
950400,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
1036800,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
1123200,32.0,5.0,3.0,5.0,2.0,0.5,0.2,210.0,100.0,10.0,2100.0,2350.0,0.0,0.0,8.1,0.1
"""


# ===========================================================================
# MAIN EXECUTION
# ===========================================================================

def main():
    """Generate all configuration files for hierarchical Mekong Delta case."""
    
    print("=" * 70)
    print("Setting up Mekong Delta Full Case (HIERARCHICAL TOPOLOGY)")
    print("=" * 70)
    print(f"Case directory: {CASE_DIR}")
    print()
    
    # Create directories
    print("1. Creating directories...")
    CASE_DIR.mkdir(parents=True, exist_ok=True)
    FORCING_DIR.mkdir(parents=True, exist_ok=True)
    print(f"   ✓ {CASE_DIR}")
    print(f"   ✓ {FORCING_DIR}")
    
    # Generate topology.csv
    print("\n2. Generating topology.csv (HIERARCHICAL)...")
    topology_content = generate_topology_csv()
    with open(CASE_DIR / "topology.csv", "w") as f:
        f.write(topology_content)
    print(f"   ✓ {len(BRANCHES)} branches defined")
    print(f"   ✓ 12 nodes (2 inlet + 5 junction + 5 outlet)")
    print("   ✓ Vam Nao: PRISMATIC (600m -> 600m) for correct physics")
    
    # Generate boundary_map.csv
    print("\n3. Generating boundary_map.csv...")
    boundary_content = generate_boundary_map_csv()
    with open(CASE_DIR / "boundary_map.csv", "w") as f:
        f.write(boundary_content)
    print(f"   ✓ 2 discharge BCs (Tien, Hau)")
    print(f"   ✓ 5 tidal level BCs (North to South)")
    
    # Generate case_config.txt
    print("\n4. Generating case_config.txt...")
    config_content = generate_case_config()
    with open(CASE_DIR / "case_config.txt", "w") as f:
        f.write(config_content)
    print("   ✓ TimeStep = 60s (stable for loops)")
    print("   ✓ Duration = 10 days + 3 warmup")
    
    # Generate biogeo_params.txt
    print("\n5. Generating biogeo_params.txt (Safe Mode)...")
    biogeo_content = generate_biogeo_params()
    with open(CASE_DIR / "biogeo_params.txt", "w") as f:
        f.write(biogeo_content)
    print("   ✓ simplified_mode = 1")
    print("   ✓ ghg_passive_mode = 1")
    
    # Generate species boundary files
    print("\n6. Generating species boundary files...")
    with open(CASE_DIR / "species_river.csv", "w") as f:
        f.write(generate_species_river())
    with open(CASE_DIR / "species_ocean.csv", "w") as f:
        f.write(generate_species_ocean())
    print("   ✓ species_river.csv (upstream)")
    print("   ✓ species_ocean.csv (downstream)")
    
    # Print network summary
    print("\n" + "=" * 70)
    print("HIERARCHICAL NETWORK TOPOLOGY")
    print("=" * 70)
    print("""
    TIEN RIVER (North, 60% flow):
      [Tan Chau] ─── Tien_Up ───┬─── Tien_Mid_A ───┬─── Tien_Mid_B ───┬─── My_Tho ─── [Sea 8]
         (1)                    │                  │                  │
                                │                  │                  └─── Ham_Luong ─── [Sea 9]
                                │                  │
                                │                  └─── Co_Chien ─────────────────────── [Sea 10]
                                │
                            Vam_Nao (PRISMATIC: 600m x 600m)
                                │
    HAU RIVER (South, 40% flow):│
      [Chau Doc] ─── Hau_Up ────┴─── Hau_Mid ─────────────────────────┬─── Dinh_An ─── [Sea 11]
         (2)                                                          │
                                                                      └─── Tran_De ─── [Sea 12]

    KEY PHYSICS:
    • Vam Nao is PRISMATIC (no Savenije convergence)
    • Flow driven by head difference dH(Tien) - dH(Hau)
    • Outlets ordered North to South with phase lag
    """)
    
    print("=" * 70)
    print("Setup Complete!")
    print("=" * 70)
    print("\nNext steps:")
    print("  1. Run: python generate_synthetic_mekong.py")
    print("  2. Build: scripts/build.bat")
    print("  3. Run: bin/Debug/CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt")
    print("  4. Audit: python audit_scientific_results.py")
    
    return 0


if __name__ == "__main__":
    exit(main())
