#!/usr/bin/env python3
"""
Setup Mekong Delta Full Case - Complete Configuration Generator
===============================================================

This script generates ALL configuration files for the complete Mekong Delta
network simulation with academically correct topology.

MEKONG DELTA TOPOLOGY (Correct Academic Representation):
=========================================================

The Mekong River splits into two main branches at Phnom Penh:
- TIEN RIVER (Tiền Giang) - Northern branch, ~60% of total flow
- HAU RIVER (Hậu Giang) - Southern branch, ~40% of total flow

These are connected by VAM NAO CANAL which redistributes flow.

                    [Tan Chau]                    [Chau Doc]
                       (1)                           (2)
                        |                             |
                   Tien_Up                        Hau_Up  
                    50 km                          45 km
                        |                             |
        [Vam Nao Jct]---(3)======VAM NAO======(4)---[Vam Nao Jct]
             Tien side   |     7 km, 18m deep  |    Hau side
                         |                     |
                    Tien_Mid                 Hau_Mid
                     60 km                    55 km
                         |                     |
                    [My Thuan]             [Can Tho]
                       (5)                    (6)
                      / | \                  /   \
                     /  |  \                /     \
              My_Tho Ham_L Co_Ch      Dinh_An   Tran_De
               70km  75km  85km        65km      60km
                 |     |     |           |          |
               (7)   (8)   (9)         (10)       (11)
            [Sea]  [Sea]  [Sea]       [Sea]      [Sea]
            
Cua      Cua Tieu  Cua   Cua Co   Cua Dinh An  Cua Tran De
         /Cua Dai  Ham   Chien
                   Luong

NODE NUMBERING (Contiguous 1-11):
- 1: Tan Chau (Tien Inlet) - Discharge BC
- 2: Chau Doc (Hau Inlet) - Discharge BC  
- 3: Vam Nao Junction - Tien side
- 4: Vam Nao Junction - Hau side
- 5: My Thuan Split (Tien distributaries)
- 6: Can Tho Split (Hau distributaries)
- 7: My Tho Mouth (Cua Tieu/Dai)
- 8: Ham Luong Mouth
- 9: Co Chien Mouth
- 10: Dinh An Mouth
- 11: Tran De Mouth

BRANCH NUMBERING (1-10):
- 1: Tien_Up (1→3) - Upper Tien, 50km
- 2: Hau_Up (2→4) - Upper Hau, 45km
- 3: Vam_Nao (3→4) - Cross-connection, 7km, DEEP (18m)
- 4: Tien_Mid (3→5) - Middle Tien, 60km
- 5: Hau_Mid (4→6) - Middle Hau, 55km
- 6: My_Tho (5→7) - Northernmost Tien distributary, 70km
- 7: Ham_Luong (5→8) - Central Tien distributary, 75km
- 8: Co_Chien (5→9) - Southernmost Tien distributary, 85km
- 9: Dinh_An (6→10) - Main Hau distributary, 65km
- 10: Tran_De (6→11) - Secondary Hau distributary, 60km

References:
- Nguyen et al. (2006) - Geometry data
- Vo Quoc Thanh (2021) - Updated bathymetry
- Fujii et al. (2003) - Vam Nao hydraulics
- MRC (Mekong River Commission) - Discharge data

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

# Project root (relative to script location)
SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent
CASE_DIR = PROJECT_ROOT / "INPUT" / "Cases" / CASE_NAME
FORCING_DIR = CASE_DIR / "forcing_data"


# ===========================================================================
# TOPOLOGY DATA (Academically Correct)
# ===========================================================================

# Node definitions
NODES = {
    1: {"name": "Tan_Chau", "type": "DISCHARGE", "description": "Tien River inlet at Vietnamese border"},
    2: {"name": "Chau_Doc", "type": "DISCHARGE", "description": "Hau River inlet at Vietnamese border"},
    3: {"name": "VamNao_Tien", "type": "JUNCTION", "description": "Vam Nao junction - Tien side"},
    4: {"name": "VamNao_Hau", "type": "JUNCTION", "description": "Vam Nao junction - Hau side"},
    5: {"name": "My_Thuan", "type": "JUNCTION", "description": "Tien distributary split point"},
    6: {"name": "Can_Tho", "type": "JUNCTION", "description": "Hau distributary split point"},
    7: {"name": "MyTho_Mouth", "type": "LEVEL", "description": "My Tho outlet (Cua Tieu/Dai)"},
    8: {"name": "HamLuong_Mouth", "type": "LEVEL", "description": "Ham Luong outlet"},
    9: {"name": "CoChien_Mouth", "type": "LEVEL", "description": "Co Chien outlet"},
    10: {"name": "DinhAn_Mouth", "type": "LEVEL", "description": "Dinh An outlet (main Hau)"},
    11: {"name": "TranDe_Mouth", "type": "LEVEL", "description": "Tran De outlet"},
}

# Branch definitions with realistic geometry
# Format: ID, Name, NodeUp, NodeDown, Length_m, Width_Up_m, Width_Down_m, Depth_m, Chezy, Group, RS
BRANCHES = [
    # ID  Name        Up Down  Length  W_Up  W_Down Depth Chezy Group RS
    (1,  "Tien_Up",    1,  3,  50000,   900,  1200,  14.0,  65,   1,  1.0),
    (2,  "Hau_Up",     2,  4,  45000,   850,  1100,  13.0,  65,   1,  1.0),
    (3,  "Vam_Nao",    3,  4,   7000,   600,   600,  18.0,  70,   3,  1.0),  # Deep canal!
    (4,  "Tien_Mid",   3,  5,  60000,  1200,  1500,  13.0,  65,   1,  1.0),
    (5,  "Hau_Mid",    4,  6,  55000,  1100,  1400,  14.0,  65,   1,  1.0),
    (6,  "My_Tho",     5,  7,  70000,  1000,  2500,  11.0,  60,   2,  1.2),
    (7,  "Ham_Luong",  5,  8,  75000,   900,  2200,  10.0,  60,   2,  1.2),
    (8,  "Co_Chien",   5,  9,  85000,   950,  2800,  12.0,  60,   2,  1.2),
    (9,  "Dinh_An",    6, 10,  65000,  1200,  3000,  11.0,  58,   2,  1.3),
    (10, "Tran_De",    6, 11,  60000,  1000,  2000,   9.0,  58,   2,  1.3),
]

# Boundary conditions
DISCHARGE_NODES = [1, 2]  # Upstream discharge BCs
LEVEL_NODES = [7, 8, 9, 10, 11]  # Downstream tidal BCs

# Forcing file names
DISCHARGE_FILES = {
    1: "Tien_Inlet.csv",
    2: "Hau_Inlet.csv",
}

TIDE_FILES = {
    7: "MyTho_Tide.csv",
    8: "HamLuong_Tide.csv",
    9: "CoChien_Tide.csv",
    10: "DinhAn_Tide.csv",
    11: "TranDe_Tide.csv",
}


# ===========================================================================
# FILE GENERATORS
# ===========================================================================

def generate_topology_csv() -> str:
    """Generate topology.csv with correct syntax."""
    lines = [
        "# ============================================================================",
        "# TOPOLOGY: Mekong Delta Full Network",
        "# ============================================================================",
        "# Academic representation of complete Vietnamese Mekong Delta",
        "#",
        "# References:",
        "#   - Nguyen et al. (2006) - Channel geometry",
        "#   - Vo Quoc Thanh (2021) - Updated bathymetry",
        "#   - Fujii et al. (2003) - Vam Nao hydraulics",
        "#   - MRC Hydrological Database - Flow splits",
        "#",
        "# Node Convention (contiguous 1-11):",
        "#   1 = Tan Chau (Tien Inlet) - Discharge BC",
        "#   2 = Chau Doc (Hau Inlet) - Discharge BC",
        "#   3 = Vam Nao Junction (Tien side)",
        "#   4 = Vam Nao Junction (Hau side)",
        "#   5 = My Thuan Split (Tien distributaries)",
        "#   6 = Can Tho Split (Hau distributaries)",
        "#   7-11 = Ocean outlets (tidal level BCs)",
        "#",
        "# Format: ID, Name, NodeUp, NodeDown, Length_m, Width_Up, Width_Down, Depth, Chezy, Group, RS",
        "# ============================================================================",
        "",
    ]
    
    # Add branch data
    for b in BRANCHES:
        bid, name, node_up, node_down, length, w_up, w_down, depth, chezy, group, rs = b
        lines.append(f"{bid}, {name}, {node_up}, {node_down}, {length}, {w_up}, {w_down}, {depth}, {chezy}, {group}, {rs}")
    
    return "\n".join(lines)


def generate_boundary_map_csv() -> str:
    """Generate boundary_map.csv with correct syntax."""
    lines = [
        "# ============================================================================",
        "# BOUNDARY CONDITIONS: Mekong Delta Full",
        "# ============================================================================",
        "# Format: Node_ID, Type, FilePath [, SpeciesName]",
        "#",
        "# Types: DISCHARGE, LEVEL, SPECIES",
        "# For SPECIES with ALL: loads all species columns from CSV",
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
        lines.append(f"# Node {node_id}: {NODES[node_id]['name']}")
        lines.append(f"{node_id}, DISCHARGE, forcing_data/{filename}")
        lines.append(f"{node_id}, SPECIES, species_river.csv, ALL")
        lines.append("")
    
    lines.append("# ===========================================================================")
    lines.append("# DOWNSTREAM TIDAL BOUNDARIES (5 RIVER MOUTHS)")
    lines.append("# ===========================================================================")
    lines.append("")
    
    # Tidal BCs
    for node_id in LEVEL_NODES:
        filename = TIDE_FILES[node_id]
        lines.append(f"# Node {node_id}: {NODES[node_id]['name']}")
        lines.append(f"{node_id}, LEVEL, forcing_data/{filename}")
        lines.append(f"{node_id}, SPECIES, species_ocean.csv, ALL")
        lines.append("")
    
    return "\n".join(lines)


def generate_case_config() -> str:
    """Generate case_config.txt."""
    config = f"""# ============================================================================
# C-GEM Network Case Configuration: {CASE_NAME}
# ============================================================================
# Complete Mekong Delta simulation with:
# - Tien River (3 distributaries: My Tho, Ham Luong, Co Chien)
# - Hau River (2 distributaries: Dinh An, Tran De)
# - Vam Nao cross-connection (bidirectional tidal flow)
#
# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
# ============================================================================

CaseName = {CASE_NAME}

# ---------------------------------------------------------------------------
# FILE PATHS (relative to this config file)
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
# CFL: dt < dx / (sqrt(g*h) + |u|)
# For dx=2000m, h=15m: dt < 2000/(12+1) ~ 150s
TimeStep = 150.0

# Grid spacing [meters]
DELXI = 2000.0

# ---------------------------------------------------------------------------
# OUTPUT OPTIONS
# ---------------------------------------------------------------------------
WriteNetCDF = 1
WriteCSV = 1
WriteReactionRates = 0
"""
    return config


def generate_biogeo_params() -> str:
    """Generate biogeo_params.txt with safe/simplified mode."""
    params = """# ============================================================================
# C-GEM Biogeochemistry Parameters: Mekong Delta (Safe Mode)
# ============================================================================
# Simplified configuration for stable simulation
# Uses 80/20 approach - bulk parameters instead of detailed kinetics
#
# SAFETY FLAGS:
#   simplified_mode = 1  -> Use bulk decay rates
#   ghg_passive_mode = 1 -> GHG as diagnostic only (no feedback)
# ============================================================================

# ---------------------------------------------------------------------------
# ENVIRONMENTAL CONDITIONS
# ---------------------------------------------------------------------------
water_temp = 28.0         # Mean temperature [°C] (tropical)
ws = 0.0005               # Settling velocity [m/s]

# ---------------------------------------------------------------------------
# LIGHT PARAMETERS
# ---------------------------------------------------------------------------
I0 = 350.0                # Surface irradiance [W/m²]
kd1 = 0.20                # Base light attenuation [1/m]
kd2_spm = 0.018           # SPM contribution [m²/mg]
kd2_phy1 = 0.008          # PHY1 self-shading
kd2_phy2 = 0.006          # PHY2 self-shading

# ---------------------------------------------------------------------------
# PHYTOPLANKTON PARAMETERS
# ---------------------------------------------------------------------------
alpha1 = 0.02             # P-I curve slope [µg C/µg Chl/W/m²]
pbmax1 = 2.5              # Max photosynthesis [1/day]
kexc1 = 0.10              # Excretion fraction
kgrowth1 = 0.10           # Growth respiration
kmaint1 = 0.025           # Maintenance [1/day]
kmort1 = 0.08             # Mortality [1/day]

alpha2 = 0.015
pbmax2 = 2.0
kexc2 = 0.10
kgrowth2 = 0.08
kmaint2 = 0.02
kmort2 = 0.06

# ---------------------------------------------------------------------------
# NUTRIENT HALF-SATURATION CONSTANTS
# ---------------------------------------------------------------------------
kdsi1 = 3.0               # DSi [µmol/L]
kn1 = 1.5                 # N [µmol/L]
kpo41 = 0.3               # PO4 [µmol/L]
kn2 = 2.0
kpo42 = 0.5

# ---------------------------------------------------------------------------
# DECOMPOSITION RATES
# ---------------------------------------------------------------------------
kox = 0.15                # Aerobic TOC decay [1/day]
kdenit = 0.04             # Denitrification [1/day]
knit = 0.12               # Nitrification [1/day]

# ---------------------------------------------------------------------------
# HALF-SATURATION CONSTANTS
# ---------------------------------------------------------------------------
ktox = 40.0               # TOC [µmol/L]
ko2 = 8.0                 # O2 for respiration [µmol/L]
ko2_nit = 4.0             # O2 for nitrification [µmol/L]
kno3 = 4.0                # NO3 for denitrification [µmol/L]
knh4 = 2.0                # NH4 for nitrification [µmol/L]
kino2 = 6.0               # O2 inhibition for denitrification [µmol/L]

# ---------------------------------------------------------------------------
# STOICHIOMETRY (Redfield)
# ---------------------------------------------------------------------------
redn = 0.151              # N:C ratio
redp = 0.0094             # P:C ratio
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
benthic_resp_20C = 70.0   # Benthic CO2 [mmol C/m²/day]
benthic_Q10 = 2.0
benthic_NH4_flux = 2.5    # Benthic NH4 [mmol N/m²/day]
benthic_PO4_flux = 0.15   # Benthic PO4 [mmol P/m²/day]

# ---------------------------------------------------------------------------
# GHG PARAMETERS (Passive Mode)
# ---------------------------------------------------------------------------
N2O_yield_nit = 0.005
N2O_yield_denit = 0.015
benthic_CH4_flux = 150.0  # [µmol/m²/day]
benthic_N2O_flux = 8.0    # [nmol/m²/day]

# ============================================================================
# SAFETY MODE FLAGS (CRITICAL!)
# ============================================================================
simplified_mode = 1       # 1 = Use simplified bulk kinetics
ghg_passive_mode = 1      # 1 = GHG diagnostic only (SAFE)
skip_bacteria = 1         # 1 = Skip explicit bacteria
skip_multipool_oc = 1     # 1 = Use single TOC pool
skip_ghg_dynamics = 0     # 0 = Calculate GHG (diagnostic)
skip_p_adsorption = 1     # 1 = Skip PO4-PIP equilibrium
"""
    return params


def generate_species_river() -> str:
    """Generate species_river.csv for upstream boundary."""
    # Header with all species
    species = """# Upstream River Water Quality: Mekong at Vietnamese Border
# =============================================================================
# Based on MRC monitoring data and field campaigns
# Format: species_index, concentration, unit
# =============================================================================

# Index, Value, Unit, Description
0, 0.1, PSU, Salinity (freshwater)
1, 15.0, ugC/L, PHY1 - Diatoms
2, 8.0, ugC/L, PHY2 - Green algae
3, 80.0, umol/L, DSi - Dissolved silica
4, 35.0, umolN/L, NO3 - Nitrate
5, 12.0, umolN/L, NH4 - Ammonium
6, 1.5, umolP/L, PO4 - Phosphate
7, 250.0, umol/L, O2 - Dissolved oxygen
8, 300.0, umolC/L, TOC - Total organic carbon
9, 80.0, mg/L, SPM - Suspended sediment
10, 1800.0, umol/L, DIC - Dissolved inorganic carbon
11, 1900.0, ueq/L, AT - Total alkalinity
12, 0.0, uatm, pCO2 - computed
13, 0.0, umol/L, CO2 - computed
14, 7.5, pH, pH
15, 0.5, umol/L, HS - Hydrogen sulfide
"""
    return species


def generate_species_ocean() -> str:
    """Generate species_ocean.csv for downstream boundary."""
    species = """# Ocean Boundary Water Quality: South China Sea
# =============================================================================
# Based on coastal monitoring data
# Format: species_index, concentration, unit
# =============================================================================

# Index, Value, Unit, Description
0, 32.0, PSU, Salinity (seawater)
1, 5.0, ugC/L, PHY1 - Marine diatoms
2, 3.0, ugC/L, PHY2 - Green algae
3, 5.0, umol/L, DSi - Dissolved silica (low)
4, 2.0, umolN/L, NO3 - Nitrate (oligotrophic)
5, 0.5, umolN/L, NH4 - Ammonium
6, 0.2, umolP/L, PO4 - Phosphate
7, 210.0, umol/L, O2 - Dissolved oxygen
8, 100.0, umolC/L, TOC - Marine DOC
9, 10.0, mg/L, SPM - Clear coastal water
10, 2100.0, umol/L, DIC
11, 2350.0, ueq/L, AT - Seawater alkalinity
12, 0.0, uatm, pCO2 - computed
13, 0.0, umol/L, CO2 - computed
14, 8.1, pH, pH (seawater)
15, 0.1, umol/L, HS - Very low
"""
    return species


# ===========================================================================
# MAIN EXECUTION
# ===========================================================================

def main():
    """Generate all configuration files for Mekong Delta Full case."""
    
    print("=" * 70)
    print("Setting up Mekong Delta Full Case")
    print("=" * 70)
    print(f"Case directory: {CASE_DIR}")
    print()
    
    # Create directories
    print("1. Creating directories...")
    CASE_DIR.mkdir(parents=True, exist_ok=True)
    FORCING_DIR.mkdir(parents=True, exist_ok=True)
    print(f"   - {CASE_DIR}")
    print(f"   - {FORCING_DIR}")
    
    # Generate topology.csv
    print("\n2. Generating topology.csv...")
    topology_content = generate_topology_csv()
    topology_path = CASE_DIR / "topology.csv"
    with open(topology_path, "w") as f:
        f.write(topology_content)
    print(f"   - {topology_path}")
    print(f"   - {len(BRANCHES)} branches defined")
    print(f"   - Nodes 1-{max(NODES.keys())} (contiguous)")
    
    # Generate boundary_map.csv
    print("\n3. Generating boundary_map.csv...")
    boundary_content = generate_boundary_map_csv()
    boundary_path = CASE_DIR / "boundary_map.csv"
    with open(boundary_path, "w") as f:
        f.write(boundary_content)
    print(f"   - {boundary_path}")
    print(f"   - {len(DISCHARGE_NODES)} discharge BCs")
    print(f"   - {len(LEVEL_NODES)} tidal level BCs")
    
    # Generate case_config.txt
    print("\n4. Generating case_config.txt...")
    config_content = generate_case_config()
    config_path = CASE_DIR / "case_config.txt"
    with open(config_path, "w") as f:
        f.write(config_content)
    print(f"   - {config_path}")
    
    # Generate biogeo_params.txt
    print("\n5. Generating biogeo_params.txt...")
    biogeo_content = generate_biogeo_params()
    biogeo_path = CASE_DIR / "biogeo_params.txt"
    with open(biogeo_path, "w") as f:
        f.write(biogeo_content)
    print(f"   - {biogeo_path}")
    print("   - simplified_mode = 1 (SAFE)")
    print("   - ghg_passive_mode = 1 (SAFE)")
    
    # Generate species boundary files
    print("\n6. Generating species boundary files...")
    
    species_river_path = CASE_DIR / "species_river.csv"
    with open(species_river_path, "w") as f:
        f.write(generate_species_river())
    print(f"   - {species_river_path}")
    
    species_ocean_path = CASE_DIR / "species_ocean.csv"
    with open(species_ocean_path, "w") as f:
        f.write(generate_species_ocean())
    print(f"   - {species_ocean_path}")
    
    # Summary
    print("\n" + "=" * 70)
    print("Setup Complete!")
    print("=" * 70)
    print()
    print("Network topology:")
    print("  TIEN RIVER SYSTEM:")
    print("    Node 1 (Tan Chau) → Tien_Up → Node 3 (Vam Nao)")
    print("    Node 3 → Tien_Mid → Node 5 (My Thuan)")
    print("    Node 5 → My_Tho → Node 7 (Sea)")
    print("    Node 5 → Ham_Luong → Node 8 (Sea)")
    print("    Node 5 → Co_Chien → Node 9 (Sea)")
    print()
    print("  HAU RIVER SYSTEM:")
    print("    Node 2 (Chau Doc) → Hau_Up → Node 4 (Vam Nao)")
    print("    Node 4 → Hau_Mid → Node 6 (Can Tho)")
    print("    Node 6 → Dinh_An → Node 10 (Sea)")
    print("    Node 6 → Tran_De → Node 11 (Sea)")
    print()
    print("  CROSS-CONNECTION:")
    print("    Node 3 ←→ Vam_Nao ←→ Node 4")
    print("    (Bidirectional tidal flow, 18m deep)")
    print()
    print("Next step: Run generate_synthetic_mekong.py to create forcing data")
    
    return 0


if __name__ == "__main__":
    exit(main())
