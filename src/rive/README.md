# C-RIVE Biogeochemistry Module

## Overview

This directory contains the **C-RIVE** (C-language RIVE) biogeochemistry module,
a comprehensive process-based model for carbon, nutrient, and oxygen cycling
in aquatic systems.

**C-RIVE is the authoritative ANSI C implementation of the Unified RIVE v1.0**,
developed at Mines Paris by the original RIVE research team (Flipo, Wang, Vilmin, Hasanyar).

---

## Module Structure

```
rive/
├── README.md               # This file
├── rive.h                  # Master include - use this to access all RIVE
├── rive_common.h           # Shared constants, temperature functions, macros
│
├── ── PARAMETERS ──────────────────────────────────────────────
├── rive_params.h           # BiogeoParams struct definition
├── rive_params.c           # Parameter loading from biogeo_params.txt
│
├── ── MAIN DRIVER ─────────────────────────────────────────────
├── biogeo.c                # Biogeo_Branch() - coordinates all processes
│
├── ── PHYTOPLANKTON ───────────────────────────────────────────
├── phytoplankton.h         # API: rive_calc_phytoplankton(), rive_temp_kin()
├── phytoplankton.c         # PHY1 (diatoms), PHY2 (greens)
│                           # - Light limitation (depth-integrated)
│                           # - Nutrient limitation (N, P, Si)
│                           # - Temperature-dependent rates
│
├── ── NUTRIENT CYCLING ────────────────────────────────────────
├── nutrients.h             # API: rive_calc_nitrification(), etc.
├── nutrients.c             # - Nitrification: NH4 → NO3
│                           # - Denitrification: NO3 → N2
│                           # - Aerobic degradation of TOC
│                           # - Silica consumption by diatoms
│
├── ── OXYGEN DYNAMICS ─────────────────────────────────────────
├── oxygen.h                # API: rive_oxygen_saturation(), etc.
├── oxygen.c                # - O2 saturation (Weiss 1970)
│                           # - Piston velocity (wind + current)
│                           # - Air-water O2 exchange
│
├── ── CARBONATE CHEMISTRY ─────────────────────────────────────
├── carbonate_chem.h        # API: crive_calc_carbonate_system()
├── carbonate_chem.c        # - DIC speciation (CO2, HCO3-, CO3--)
│                           # - pH calculation (Newton-Raphson)
│                           # - pCO2 calculation
│                           # - CO2 air-water flux
│                           # - k600 parameterizations
│
├── ── GREENHOUSE GASES ────────────────────────────────────────
├── ghg_module.h            # API: crive_calc_ghg_system()
├── ghg_module.c            # - N2O: from nitrification + denitrification
│                           # - CH4: methanogenesis, oxidation, ebullition
│                           # - 2-step nitrification (NH4 → NO2 → NO3)
│
└── ── NUMERICAL SOLVER ────────────────────────────────────────
    rk4_solver.h            # API: crive_rk4_step()
    rk4_solver.c            # 4th-order Runge-Kutta with adaptive stepping
```

---

## State Variables (17 species)

| Index | Symbol | Name | Unit | Description |
|-------|--------|------|------|-------------|
| 0 | SAL | Salinity | PSU | Conservative tracer |
| 1 | PHY1 | Diatoms | µgC/L | Siliceous phytoplankton |
| 2 | PHY2 | Green algae | µgC/L | Non-siliceous phytoplankton |
| 3 | DSi | Dissolved silica | µmol/L | Limiting nutrient for diatoms |
| 4 | NO3 | Nitrate | µmol/L | Oxidized nitrogen |
| 5 | NH4 | Ammonium | µmol/L | Reduced nitrogen (preferred by algae) |
| 6 | PO4 | Phosphate | µmol/L | Often limiting in freshwater |
| 7 | O2 | Dissolved oxygen | µmol/L | ~250 µmol/L at saturation |
| 8 | TOC | Total organic carbon | µmol/L | Substrate for heterotrophs |
| 9 | SPM | Suspended matter | mg/L | Affects light, adsorbs P |
| 10 | DIC | Dissolved inorganic C | µmol/L | CO2 + HCO3- + CO3-- |
| 11 | AT | Total alkalinity | µmol/L | Acid-neutralizing capacity |
| 12 | pCO2 | Partial pressure CO2 | µatm | Diagnostic - drives flux |
| 13 | CO2 | Aqueous CO2 | µmol/L | Diagnostic |
| 14 | pH | pH | - | Diagnostic |
| 15 | HS | Hydrogen sulfide | µmol/L | Anoxic indicator |
| 16 | ALKC | Carbonate alkalinity | µmol/L | Diagnostic (iterations) |

---

## Reaction Rates (27 tracked)

| Category | Reactions |
|----------|-----------|
| **Primary Production** | GPP_1, GPP_2, NPP_NO3, NPP_NH4, NPP (total) |
| **Phytoplankton Death** | PHY_DEATH_1, PHY_DEATH_2, PHY_DEATH (total) |
| **Nutrient Cycling** | NIT (nitrification), DENIT, SI_CONS |
| **Organic Matter** | AER_DEG (aerobic degradation) |
| **Oxygen Exchange** | O2_EX, O2_EX_S (surface flux) |
| **Carbon Exchange** | CO2_EX, CO2_EX_S, DIC_REACT, TA_REACT |
| **GHG** | N2O_NIT, N2O_DENIT, N2O_EX, CH4_PROD, CH4_OX, CH4_EX, CH4_EBUL |
| **Sediment** | EROSION, DEPOSITION |

---

## Key Features

### 1. Phytoplankton Module (`phytoplankton.c`)

**Two functional groups:**
- **PHY1 (Diatoms)**: Require silica, dominate in turbulent/nutrient-rich water
- **PHY2 (Green algae)**: No silica requirement, dominate in stratified water

**Processes:**
```
Light → Gross Primary Production (GPP)
GPP - Excretion - Growth cost - Maintenance = Net Primary Production (NPP)
NPP uses either NO3 or NH4 (preference for NH4)
Mortality → Organic matter pool
```

**Key functions:**
- `rive_light_limitation()` - Depth-integrated P-I curve
- `rive_temp_kin()` - Temperature correction (Q10, Arrhenius)
- `rive_calc_phytoplankton()` - All growth calculations

### 2. Nutrient Module (`nutrients.c`)

**Nitrogen cycle:**
```
       ┌─────────────────┐
       │   Nitrification │
NH4 ──►│ (O2 required)   │──► NO3
       └─────────────────┘     │
                               ▼
       ┌─────────────────┐     
       │ Denitrification │◄────┘
NO3 ──►│ (low O2, TOC)   │──► N2 (lost)
       └─────────────────┘
```

**Key functions:**
- `rive_calc_nitrification()` - NH4 → NO3 (produces N2O)
- `rive_calc_denitrification()` - NO3 → N2 (produces N2O)
- `rive_calc_aerobic_degradation()` - TOC mineralization

### 3. Oxygen Module (`oxygen.c`)

**O2 budget:**
```
Production:  Photosynthesis (NPP)
Consumption: Respiration, nitrification, TOC degradation
Exchange:    Air-water (reaeration when undersaturated)
```

**Key functions:**
- `rive_oxygen_saturation()` - Weiss (1970) temperature/salinity correction
- `rive_piston_velocity()` - k = f(wind, current, depth)
- `rive_calc_o2_exchange()` - Flux = k/H × (O2_sat - O2)

### 4. Carbonate Chemistry (`carbonate_chem.c`)

**DIC speciation:**
```
                  K1              K2
CO2(aq) + H2O ←→ HCO3- + H+ ←→ CO3-- + 2H+
```

**pH calculation:**
- Newton-Raphson iteration on alkalinity equation
- Includes borate, sulfide contributions

**CO2 flux:**
- `k600` parameterization options: Strahler, Raymond, Abril
- Adjusts for Schmidt number at temperature

### 5. GHG Module (`ghg_module.c`)

**N2O sources:**
- Nitrification byproduct: ~0.4% of NH4 oxidized
- Incomplete denitrification: ~1% of NO3 reduced

**CH4 dynamics:**
- Methanogenesis in anoxic sediments
- Aerobic oxidation in water column
- Ebullition when CH4 > threshold

---

## How to Use

### Basic Usage

```c
#include "rive/rive.h"  // Master include

// Initialize parameters (from file or defaults)
LoadBiogeoParams("INPUT/Cases/MyCase/biogeo_params.txt");
InitializeBiogeoParameters(branch);

// In time loop:
Biogeo_Branch(branch, dt);      // Main biogeochemistry
Biogeo_GHG_Branch(branch, dt);  // Optional: GHG dynamics
```

### Adding a New Process

1. **Choose the right file** based on process type:
   - Phyto-related → `phytoplankton.c`
   - Nutrient-related → `nutrients.c`
   - Gas exchange → `oxygen.c` or `ghg_module.c`

2. **Add function to header** with documentation:
   ```c
   /**
    * Calculate new process rate
    * @param conc Concentration [µmol/L]
    * @param temp Temperature [°C]
    * @return Rate [µmol/L/s]
    */
   double rive_calc_new_process(double conc, double temp);
   ```

3. **Implement in .c file** with literature reference:
   ```c
   /* Reference: Author et al. (Year) Journal */
   double rive_calc_new_process(double conc, double temp) {
       double rate = ...;
       return rate;
   }
   ```

4. **Call from biogeo.c** in `Biogeo_Branch()`:
   ```c
   double new_rate = rive_calc_new_process(conc, temp);
   branch->reaction_rates[CGEM_REACTION_NEW][i] = new_rate;
   ```

---

## Parameter File Format

`biogeo_params.txt` uses key = value format:

```
# Phytoplankton
pbmax1 = 2.0         # Max photosynthesis PHY1 [/day]
pbmax2 = 1.5         # Max photosynthesis PHY2 [/day]
kmort1 = 0.05        # Mortality PHY1 [/day]

# Nutrients
kox = 0.1            # TOC oxidation rate [/day]
knit = 0.1           # Nitrification rate [/day]
kdenit = 0.05        # Denitrification rate [/day]

# Gas exchange
wind_speed = 3.5     # Mean wind at 10m [m/s]
pco2_atm = 420.0     # Atmospheric pCO2 [µatm]
```

---

## Citations

**Main C-RIVE reference:**
```bibtex
@article{wang2018crive,
  title={Time-dependent global sensitivity analysis of the {C-RIVE} biogeochemical model},
  author={Wang, Shuaitao and Flipo, Nicolas and Romary, Thomas},
  journal={Water Research},
  volume={144},
  pages={341--355},
  year={2018}
}
```

**Unified RIVE:**
```bibtex
@article{hasanyar2022unified,
  title={Unified {RIVE} v1.0: a multi-environment aquatic biogeochemical model},
  author={Hasanyar, Masihullah and Flipo, Nicolas and Vilmin, Lauriane and Wang, Shuaitao},
  journal={Biogeosciences},
  year={2022}
}
```

**Process-specific:**
| Process | Reference |
|---------|-----------|
| RIVE foundations | Billen et al. (1994) |
| Carbonate equilibrium | Zeebe & Wolf-Gladrow (2001) |
| CO2 solubility | Weiss (1974) |
| N2O from nitrification | Cébron et al. (2005) |
| N2O from denitrification | Garnier et al. (2007) |
| CH4 in estuaries | Borges & Abril (2011) |
